#include <iostream>
#include <algorithm> 
#include <fstream>
#include <string>
#include <vector>
#include "HolographicVQCD.h"
#include "DeepInelasticScattering.h"
#include "F2.h"
#include "Spectra.h"
#include "Kernel.h"
#include "Reggeon.h"
#include "U1NNMode.h"
#include "methods/search.hpp"
#include "methods/interpolation/Poly_Interp.hpp"
#include "methods/vectorOperators.hpp"

class F2IzNIntegrand: public IzNIntegrand
{
    public:
        F2IzNIntegrand(const Poly_Interp<double> &f1, const Poly_Interp<double> &f2, const U1NNMode &f3, const Poly_Interp<double> &f4);
        double operator()(const double x);
};

F2IzNIntegrand::F2IzNIntegrand(const Poly_Interp<double> &f1, const Poly_Interp<double> &f2, const U1NNMode &f3, const Poly_Interp<double> &f4):
                IzNIntegrand(f1, f2, f3, f4) {}

double F2IzNIntegrand::operator()(const double x) {return func1.interp(x) * func2.interp(x) * func3.factor(x) * func4.interp(x);}


extern"C"
{
    void dqags_(double (*f)(double *, void *), void * params,
                double * a, double * b, double * epsabs, double * epsrel, 
                double * result, double * abserr, int * neval,int * ier,
                int * limit, int * lenw, int * last, int * iwork, double * work);
}

F2::F2(std::string file_path) : DeepInelasticScattering(file_path)
{
    std::cout << "Loaded F2." << std::endl;
}

double F2::IzN(const std::vector<double> &kin, const Reggeon &reg)
{
    // Get the kinematical values and J
    const double Q2 = kin[0];
    const double J = reg.getJ();
    const std::string reg_name = reg.getName();
    std::vector<std::vector<double> > wf = reg.getWf();
    // Define e^(As(2 - jn))
    std::vector<double> fact1 = exp((2.0-J) * Astring);
    Poly_Interp<double> f1(z, fact1, 4);
    // Search for the correct mode
    U1NNMode mode = searchMode(Q2);
    Poly_Interp<double> f4(wf[0], wf[1], 4);
    // Choose the background potentials factor accordingly
    Poly_Interp<double> bckPotFac;
    if(reg_name == "gluon") bckPotFac = GluonPotFactor;
    else bckPotFac = MesonPotFactor;
    // Define the integrand object
    F2IzNIntegrand integrand(f1, bckPotFac, mode, f4);
    void * params = &integrand;
    // Compute the integral
    double izn = 0.0, abserr = 0.0;
    // If we look at the potential profiles vs u we see that this is good enough
    double a = z[0], b = z.back();
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