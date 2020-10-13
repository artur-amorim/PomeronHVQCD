#include <iostream>
#include <cmath>
#include <fstream>
#include "HolographicVQCD.h"
#include "U1NNMode.h"
#include "methods/interpolation/Poly_Interp.hpp"
#include "methods/vectorOperators.hpp"


// Fortran functions needed to compute the Nonnormalizable modes
extern"C"
{
    void colnew_(int * NCOMP, int * M, double * ALEFT, double * ARIGHT, double * Zeta, int * IPAR, int * LTOL,
                double * TOL, double * FIXPNT, double * ISPACE, double * FSPACE, int * IFLAG,
                void (*f) (double *, double *, double *, double *),
                void (*df) (double *, double *, double *, double *),
                void (*g) (int *, double *, double *, double *),
                void (*dg) (int *, double *, double *, double *),
                void (*guess) (double *, double *, double *), double * PARS);

    void appsln_ (double * X, double * Z, double * FSPACE, double * ISPACE);
}

std::vector<double> U1NNMode::u = {};
Poly_Interp<double> U1NNMode::t1 = Poly_Interp<double>({},{},4);

void setupU1NNcomputation()
{
    // Get u
    U1NNMode::u = hvqcdU1NNMode().u();
    const int n = U1NNMode::u.size();
    // t1Y = -( 1 + 2 dwdPhi dPhidA / w + dtaudA Vf(0,1)/Vf + dPhidA Vf(1,0)/Vf ) exp(A) / (G q)
    std::vector<double> A = hvqcdU1NNMode().A();
    std::vector<double> q = hvqcdU1NNMode().q();
    std::vector<double> Phis = hvqcdU1NNMode().Phi();
    std::vector<double> taus = hvqcdU1NNMode().tau();
    std::vector<double> dPhidA = hvqcdU1NNMode().dPhi();
    std::vector<double> dtaudA = hvqcdU1NNMode().dtaudA();
    std::vector<double> G = hvqcdU1NNMode().G();
    double a1 = hvqcdU1NNMode().get_a1(), a2 = hvqcdU1NNMode().get_a2();

    double dwdPhi, w, dlogVfdPhi, dlogVfdtau;
    std::vector<double> t1Y(n);
    for(int i = 0; i < n; i++)
    {
        w = hvqcdU1NNMode().w(Phis[i]); dwdPhi = hvqcdU1NNMode().dwdPhi(Phis[i]);
        dlogVfdPhi =  hvqcdU1NNMode().dVf0dPhi(Phis[i]) / hvqcdU1NNMode().Vf0(Phis[i]);
        dlogVfdtau =  -2*taus[i]*(-a1 + a2 + a1*a2*taus[i]*taus[i])/(1+a1*taus[i]*taus[i]);
        // Compute t1Y
        t1Y[i] = -(1+2*dwdPhi * dPhidA[i] / w + dtaudA[i] * dlogVfdtau + dPhidA[i] *  dlogVfdPhi) * std::exp(A[i]) / (G[i] * q[i]);
    }

    std::reverse(std::begin(U1NNMode::u), std::end(U1NNMode::u));
    std::reverse(std::begin(t1Y), std::end(t1Y));

    U1NNMode::t1 = Poly_Interp<double>(U1NNMode::u, t1Y, 4);
}

U1NNMode::U1NNMode(const double qq2):
    q2(qq2),
    ISPACE(new double[IIDIM]), FSPACE(new double[IFDIM])
    {}

U1NNMode::U1NNMode(const U1NNMode &mode): q2(mode.q2),
                                          ISPACE(nullptr), FSPACE(nullptr)
{
    ISPACE = new double[IIDIM];
    for(int i = 0; i < IIDIM; i++) ISPACE[i] = mode.ISPACE[i];
    FSPACE = new double[IFDIM];
    for(int i = 0; i < IFDIM; i++) FSPACE[i] = mode.FSPACE[i];
}

void U1NNMode::copy(const U1NNMode &rhs)
{
    q2 = rhs.q2;
    delete[] ISPACE; delete[] FSPACE;
    ISPACE = new double[IIDIM]; FSPACE = new double[IFDIM];
    for(int i = 0; i < IIDIM; i++) ISPACE[i] = rhs.ISPACE[i];
    for(int i = 0; i < IFDIM; i++) FSPACE[i] = rhs.FSPACE[i];
}

double U1NNMode::Q2() const
{
    // Return Q2
    return this->q2;
}

double U1NNMode::fQ(const double x) const
{
    /*
        Returns the value of fQ given x <= z.back(). If x > z.back() throw runtime_error
    */
    if (x > u.back()) std::runtime_error("U1NNMode::fQ: x out of range");
    double * X = new double(x);
    double * Z = new double[2];
    appsln_(X, Z, FSPACE, ISPACE);
    double fq = Z[0];
    // Delete X and Z
    delete X;
    delete[] Z;
    return fq;
}
double U1NNMode::dfQ(const double x) const
{
    /*
        Returns the value of dfQ/du given x <= u.back(). If x > u.back() throw runtime_error
    */
    if (x > u.back()) std::runtime_error("U1NNMode::dfQ: x out of range");
    double * X = new double(x);
    double * Z = new double[2];
    appsln_(X, Z, FSPACE, ISPACE);
    double dfq = Z[1];
    // Delete X and Z
    delete X;
    delete[] Z;
    return dfq;
}

double U1NNMode::factor(const double x) const
{
    /*
        Returns the value of fQ^2 + (dfQ/du)^2/Q2 given x <= u.back(). If x > u.back() throw runtime_error
    */
    if (x > u.back()) throw std::runtime_error("U1NNMode::factor: x is out of range");
    if (q2 < 1e-9) return 1;
    double * X = new double(x);
    double * Z = new double[2];
    appsln_(X, Z, FSPACE, ISPACE);
    double fq = Z[0];
    double dfq = Z[1];
    double fact = fq*fq + dfq*dfq / q2; 
    // Delete X and Z
    delete X;
    delete[] Z;
    return fact;
}

void U1NNMode::f(double * X, double * Z, double *F, double * PARS)
{
    double q2 = *PARS;
    F[0] = Z[1];
    F[1] = q2 * Z[0] + t1.interp(*X) * Z[1];
    return ;
}

void U1NNMode::df(double * X, double * Z, double * DF, double * PARS)
{
    // REMEMBER C++ Array is row ordered while Fortran is column ordered
    // This function is parsed to Fortran code so make sure DF is column ordered
    double q2 = *PARS;
    DF[0] = 0.0;
    DF[1] = q2;
    DF[2] = 1.0;
    DF[3] = t1.interp(*X); 
    return ;
}

void U1NNMode::g(int * I, double * Z, double * G, double * PARS)
{
    if(*I == 1) *G = Z[0] - 1.0;
    if(*I == 2) *G = Z[0] ;
    return ;
}

void U1NNMode::dg(int * I, double * Z, double * DG, double * PARS)
{
    /*if (*I == 1)
    {
        DG[0] = 1.0;
        DG[1] = 0.0;
    }
    else
    {
        DG[0] = 0.0;
        DG[1] = 1.0;
    }
    */
   DG[0] = 1.0;
   DG[1] = 0.0;
    return ;
}

void U1NNMode::guess(double * X, double * Z, double * DMVAL)
{
    return;
}

void U1NNMode::computeMode()
{
    // Check if everything is ready to start the computation
    if (U1NNMode::u.size() == 0) setupU1NNcomputation();
    // Create necessary variables
    int * NCOMP = new int(2);
    int * M = new int[2];
    M[0] = 1; M[1] = 1;
    double * ALEFT = new double(u[0]);
    double * ARIGHT = new double(u.back());
    double * Zeta = new double[2];
    Zeta[0] = *ALEFT; Zeta[1] = *ARIGHT;
    int * NOTOL = new int(2);
    int * IPAR = new int[13];
    IPAR[0] = 0;
    IPAR[1] = 0;
    IPAR[2] = 5;
    IPAR[3] = *NOTOL;
    IPAR[4] = IFDIM;
    IPAR[5] = IIDIM;
    IPAR[6] = 1;
    IPAR[7] = 0;
    IPAR[8] = 0;
    IPAR[9] = 0;
    IPAR[10] = 0;
    IPAR[11] = 0;
    IPAR[12] = 0;
    int * LTOL = new int[2];
    LTOL[0] = 1; LTOL[1] = 2;
    double * TOL = new double[2];
    TOL[0] = 1e-15; TOL[1] = 1e-15;
    double * FIXPNT = new double[1];
    int * IFLAG = new int();

    double * PARS = new double (q2);
    colnew_(NCOMP, M, ALEFT, ARIGHT, Zeta, IPAR, LTOL, TOL, FIXPNT, ISPACE, FSPACE, IFLAG, f, df, g, dg, guess, PARS);
    // Delete the variables
    delete NCOMP;
    delete[] M;
    delete ALEFT;
    delete ARIGHT;
    delete[] Zeta;
    delete NOTOL;
    delete[] IPAR;
    delete[] LTOL;
    delete[] TOL;
    delete[] FIXPNT;
    delete IFLAG;
    delete PARS;

    return ;
}

void U1NNMode::saveMode(std::string file_path) const
{
    /*
        Given a string with the file path it saves the NN modes quantities in a file.
    */
    // Create the file
    std::ofstream myfile;
    if(file_path == "")
    {
        std::string file_name;
        std::cout << "Please insert a file path first: ";
        std::cin >> file_name;
        myfile.open(file_name);
    }
    else
    {
        myfile.open(file_path);
    }
    // Write the value of Q2
    myfile << "Q2:" << '\t' << q2 << std::endl;
    myfile << std::endl;
    // Write the values of u, fQ, dfQ and factor in the file
    myfile << "u" << '\t' << "fQ" << '\t' << "dfQ/du" << '\t' << "fQ^2 + (dfQ/du)^2/Q2" << std::endl;
    for(int i = 0; i < u.size(); i++)
    {
        myfile << u[i] << '\t' << fQ(u[i]) << '\t' << dfQ(u[i]) << '\t' << factor(u[i]) << std::endl;
    }
    // Close the file
    myfile.close();
}

U1NNMode& U1NNMode::operator= (const U1NNMode &rhs)
{
    if (this == &rhs) return *this;
    copy(rhs);
    return *this ;
}

bool U1NNMode::operator< (const U1NNMode &mode) const
{
    return this->q2 < mode.q2;
}

bool U1NNMode::operator> (const U1NNMode &mode) const
{
    return this->q2 > mode.q2;
}

U1NNMode::~U1NNMode()
{
    delete[] ISPACE;
    delete[] FSPACE;
}