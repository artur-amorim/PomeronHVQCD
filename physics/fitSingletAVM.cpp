#include <iostream>
#include <exception>
#include "HolographicVQCD.h"
#include "schrodinger/chebspec.h"
#include "schrodinger/schrodinger.h"
#include "methods/optimization/NelderMead.hpp"

using namespace std;

struct Function
{
    HVQCD hvqcd;
    vector<double> VAVM;
    double rho_mass;
    Function(const double ssc, const double kksc, const double wwsc, const double WW0, const double ww0,
            const double kkU1, const double wwU1, const double VVgIR, const double WWIR, const double kkIR,
            const double wwIR, const double WW1, const double kk1, const double ww1, const double xxf,
            const double ttau0);
    double operator() (const vector<double> &X);
};

Function::Function(const double sc, const double ksc, const double wsc, const double W0, const double w0,
            const double kU1, const double wU1, const double VgIR, const double WIR, const double kIR,
            const double wIR, const double W1, const double k1, const double w1, const double xf,
            const double tau0): hvqcd(sc, ksc, wsc, W0, w0, kU1, wU1, VgIR, WIR, kIR, wIR, W1, k1, w1, xf, tau0, 0, 0),
            VAVM(vector<double>()), rho_mass(0)
            {
                try
                {
                    // Prepare the background
                    hvqcd.solve();
                    // Now compute VVM and VAVM
                    vector<double> VVM = computeVectorMesonPotential(hvqcd);
                    VAVM = computeAxialVectorMesonNonSingletPotential(hvqcd, VVM);
                    // Compute the mass of the rho meson
                    vector<double> us = hvqcd.u(), VMMasses = computeMasses(us, VVM, 6, "cheb");
                    rho_mass = VMMasses[0];
                }
                catch(...)
                {
                    throw(runtime_error("Can not initialize Function struct! Aborting fit."));
                }
            }

double Function::operator() (const vector<double> &X)
{
    double Za = X[0], ca = X[1];

    cout << "Za: " << Za << " ca: " << ca << endl;
    
    hvqcd.setZa(Za); hvqcd.setca(ca);

    // Compute the potentials
    vector<double> VSingletAVM;
    try
    {
        VSingletAVM = computeAxialVectorMesonSingletPotential(hvqcd, VAVM);
    }
    catch(...)
    {
        cout << "Potential not computed, check this" << endl;
        cout << "erms: " << 1e99 << endl;
        return 1e99;
    }

    // Compute the masses
    vector<double> SingletAVMMasses;
    try
    {
        vector<double> us = hvqcd.u();
        SingletAVMMasses = computeMasses(us, VSingletAVM, 2, "cheb");
        if(SingletAVMMasses.size() == 0) throw(runtime_error("Only negative values for Singlet AVM masses"));
    }
    catch(...)
    {
        cout << "Unable to compute spectrum. Check these values" << endl;
        cout << "erms: " << 1e99 << endl;
        return 1e99;
    }
    double erms = 0;
    // In case we have negative m2 the mass container will be empty. In that case we fill it with zero entries
    if( SingletAVMMasses.size() == 0) SingletAVMMasses = vector<double>(2,0);
    for(int i = 0; i < Rf1_rho.size(); i++) erms += fabs((SingletAVMMasses[i]/rho_mass-Rf1_rho[i])/Rf1_rho[i]);
   
    if (std::isnan(erms))
    {
        erms = 1e99;
        cout << "None erms" << endl;
        cout << "erms: " << erms << endl;
        return erms;
    }
    cout << erms << endl;
    return erms;
}

int main(int argc, char ** argv)
{
    double sc, ksc, wsc, W0, w0, kU1, wU1;
    double VgIR, WIR, kIR, wIR, W1, k1, w1;
    double xf = 2.0/3, tau0, Za, ca;
    if (argc < 16)
    {
        sc = 3.0; ksc = 3.0; wsc = 1.56; W0 = 2.5; w0 = 1.26; kU1 = 11./9; wU1 = 0.0;
        VgIR = 2.05; WIR = 0.9; kIR = 1.8; wIR = 5.0; W1 = 0.0; k1 = -0.23;
        w1 = 0.0; tau0 = 1.; Za = 1; ca = 1;
    }
    else
    {
        sc = stod(argv[1]); ksc = stod(argv[2]); wsc = stod(argv[3]); W0 = stod(argv[4]); w0 = stod(argv[5]);
        kU1 = stod(argv[6]);
        wU1 = stod(argv[7]); VgIR = stod(argv[8]); WIR = stod(argv[9]); kIR = stod(argv[10]); wIR = stod(argv[11]);
        W1 = stod(argv[12]); k1 = stod(argv[13]); w1 = stod(argv[14]); tau0 = stod(argv[15]);
        Za = stod(argv[16]); ca = stod(argv[17]);
    }
    
    cout << "Starting fit with values" << endl;
    cout << "Za: " << Za << " ca: " << ca << endl;

    // Fit the model to the spectrum
    chebSetN(800);

    vector<double> x_guess = {Za, ca};
    vector<double> deltas = {0.1, 0.1};

    Function func(sc, ksc, wsc, W0, w0, kU1, wU1, VgIR, WIR, kIR, wIR, W1, k1, w1, xf, tau0);
    
    vector<double> xop = optimFunction(x_guess, func, deltas);
    
    // Show the optimum values found
    Za = xop[0]; ca = xop[1];

    HVQCD hvqcd(sc, ksc, wsc, W0, w0, kU1, wU1, VgIR, WIR, kIR, wIR, W1, k1, w1, xf, tau0, Za, ca);
    hvqcd.solve();

    // Computing the mass ratios
    computeHVQCDRatios(hvqcd);

    double erms = func(xop);
    cout << "Minimum erms found for ";
    cout << "Za: " << Za << " ca: " << ca << " erms: " << erms << endl;
    cout << "Parameters used for the background" << endl;
    cout << "sc: " << sc << " ksc: " << ksc << " wsc: " << wsc << " W0: " << W0 << " w0: " << w0 << " kU1: " << kU1;
    cout << " wU1: " << wU1 << " VgIR: " << VgIR << " WIR: " << WIR << " kIR: " << kIR << " wIR: " << wIR << " W1: " << W1;
    cout << " k1: " << k1 << " w1: " << w1 << " tau0: " << tau0 << endl;

    return 0;
}