#include <iostream>
#include <cmath>
#include <boost/math/special_functions/bessel.hpp>

using namespace std;
using namespace boost::math;

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

void f(double * X, double * Z, double *F, double * PARS)
{
    double q2 = *PARS;
    F[0] = Z[1];
    F[1] = q2 * Z[0] + Z[1] / (*X);
    return ;
}

void df(double * X, double * Z, double * DF, double * PARS)
{
    // REMEMBER C++ Array is row ordered while Fortran is column ordered
    // This function is parsed to Fortran code so make sure DF is column ordered
    double q2 = *PARS;
    DF[0] = 0.0;
    DF[1] = q2;
    DF[2] = 1.0;
    DF[3] = 1 / (*X); 
    return ;
}

void g(int * I, double * Z, double * G, double * PARS)
{
    if(*I == 1) *G = Z[0] - 1.0;
    if(*I == 2) *G = Z[0] ;
    return ;
}

void dg(int * I, double * Z, double * DG, double * PARS)
{
    if (*I == 1)
    {
        DG[0] = 1.0;
        DG[1] = 0.0;
    }
    else
    {
        DG[0] = 0.0;
        DG[1] = 1.0;
    }
    return ;
}

void guess(double * X, double * Z, double * DMVAL)
{
    return;
}

int main(int argc, char ** argv)
{
    double Q2 = stod(argv[1]);
    const double a = 1e-5;
    const double b = 24.8733 / sqrt(Q2);
    // Create necessary variables
    const int IFDIM = 10000;
    const  int IIDIM = 10000;
    double * ISPACE = new double[IIDIM];
    double * FSPACE = new double[IFDIM];
    int * NCOMP = new int(2);
    int * M = new int[2];
    M[0] = 1; M[1] = 1;
    double * ALEFT = new double(1e-5);
    double * ARIGHT = new double(2*b);
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

    double * PARS = new double (Q2);
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

    // Now that the mode is computed let's compare the numerical solution with the analytic one

    delete[] ISPACE;
    delete[] FSPACE;

    for(double x = a; x <= b; x += (b-a)/200)
    {
        double * X = new double(x);
        double * Z = new double[2];
        appsln_(X, Z, FSPACE, ISPACE);
        double fq = Z[0];
        double exact_sol = x * sqrt(Q2) * cyl_bessel_k(1, x * sqrt(Q2));
        // Delete X and Z
        delete X;
        delete[] Z;
        cout << x << '\t' << fq << '\t' << exact_sol << '\t' << fabs(fq - exact_sol) / exact_sol << endl;
    }

    return 0;
}