#include "HolographicVQCD.h"
#include "FL.h"
#include "methods/interpolation/Poly_Interp.hpp"
#include "methods/vectorOperators.hpp"

FLIzNIntegrand::FLIzNIntegrand(const Poly_Interp<double> &f1, const Poly_Interp<double> &f2, const U1NNMode &f3, const Poly_Interp<double> &f4):
                IzNIntegrand(f1, f2, f3, f4) {}

double FLIzNIntegrand::operator()(const double x)
{
    double dfq = func3.dfQ(x);
    return func1.interp(x) * func2.interp(x) * (dfq * dfq / func3.Q2()) * func4.interp(x);
}

extern"C"
{
    void dqags_(double (*f)(double *, void *),  void * params,
                double * a, double * b, double * epsabs, double * epsrel, 
                double * result, double * abserr, int * neval,int * ier,
                int * limit, int * lenw, int * last, int * iwork, double * work);
}


FL::FL(const bool rrsslog, const std::string file_path) : DeepInelasticScattering(rrsslog, file_path)
{
    std::cout << "Loaded FL." << std::endl;
}

double FL::IzN(const std::vector<double> &kin, const Reggeon &reg)
{
    // Get the kinematical values and J
    const double Q2 = kin[0];
    const double J = reg.getJ();
    const std::string reg_name = reg.getName();
    std::vector<std::vector<double> > wf = reg.getWf();
    // Define e^(As(1.5 - jn)) or e^(As(2-j_n))
    std::vector<double> fact1;
    if(reg_name == "gluon") fact1 = exp((1.5-J) * Astring);
    else fact1 = exp((2.0-J) * Astring);
    Poly_Interp<double> f1(u, fact1, 4);
    // Search for the correct mode
    U1NNMode mode = searchMode(Q2);
    // Interpolate the wavefunction after computing the respective u value.
    if(reg_name == "gluon") for(int i = 0; i < wf[0].size(); i++) wf[0][i] = ufunc.interp(wf[0][i]);
    Poly_Interp<double> f4(wf[0], wf[1], 4);
    // Choose the background potentials factor accordingly
    Poly_Interp<double> bckPotFac;
    if(reg_name == "gluon") bckPotFac = potFactor;
    else bckPotFac = MesonPotFactor;
    // Define the integrand object
    FLIzNIntegrand integrand(f1, bckPotFac, mode, f4);
    void * params = &integrand;
    // Compute the integral
    double izn = 0.0, abserr = 0.0;
    double a = u[0], b = u.back();
    // Absolute and relative tolerances desidered
    double epsabs = 1e-9, epsrel = 1e-9;
    // Output variables 
    int neval = 0, ier;
    // Setup the workspace for the method
    int limit = 100000;
    int lenw = 1000000;
    int last = 0;
    int * iwork = new int[limit];
    double * work = new double[lenw];
    // Evaluate the integral
    dqags_(f, params, &a, &b, &epsabs, &epsrel, &izn, &abserr, &neval, &ier, &limit, &lenw, &last, iwork, work);
    izn = std::pow(Q2, J) * izn;
    // Free workspace from memory
    delete[] iwork;
    delete[] work;
    return izn;
}