#include <iostream>
#include <string>
#include <vector>
#include "HolographicVQCD.h"
#include "HQCDP.h"
#include "GluonKernel.h"
#include "MesonKernel.h"
#include "SigmaProtonProton.h"
#include "methods/optimization/NelderMead.hpp"
#include "schrodinger/chebspec.h"

using namespace std;

int main(int argc, char ** argv)
{
    string data_path;
    double invls, g1, g2, m1;
    if(argc < 6)
    {
        data_path = "expdata/SigmaProtonProton/SigmaProtonProton_data_W_lt_10000.txt";
        invls = 0.25; g1 = 0; g2 = 0; m1 = 0;
    }
    else
    {
        data_path = argv[1]; invls = stod(argv[2]); g1 = stod(argv[3]); g2 = stod(argv[4]); m1 = stod(argv[5]);
    }
    cout << "Starting the fit with" << endl;
    cout << "invls: " << invls << " g1: " << g1 << " g2: " << g2 << " m1: " << m1 << endl;

    double mq = hvqcd().QuarkMass();
    SigmaProtonProton sigma_pp(data_path);

    // Setup Gluon and Meson Kernels and GNs vector
    vector<double> gluon_pars = {invls, 0, 0, 0, 0, 0, 0};
    GluonKernel gluon(2, gluon_pars);
    vector<double> meson_pars = gluon_pars;
    MesonKernel meson(1, meson_pars);
    vector<double> GNs = {g1, g2, m1};

    // Setup HQCDP object
    HQCDP hqcdp;
    hqcdp.addProcessObservable(sigma_pp);
    hqcdp.addKernel(gluon);
    hqcdp.addKernel(meson);
    hqcdp.setGNs(GNs);

    // initialise Chebyschev matrices
    chebSetN(1000);

    vector<double> x_guess = {invls, g1, g2, m1};

       // Definition of the function we want to fit
    auto func = [&hqcdp] (const vector<double> &x)
    {
        // gluon pars correspond to x(0)
        vector<double> gluon_pars = {x[0], 0, 0, 0, 0, 0, 0};
        // meson pars correspond are the same as gluon_pars
        vector<double> meson_pars = gluon_pars;
        vector<vector<double> > kernels_pars = {gluon_pars, meson_pars};
        // GNs correspond to x(1), x(2), x(3)
        vector<double> gpars = {x[1], x[2], x[3]};
        // Update the spectrum
        hqcdp.computeSpectrum(kernels_pars);
        // Update the GNs
        hqcdp.setGNs(gpars);
        // Compute chi2
        double CHI2 = hqcdp.chi2();
        if(std::isnan(CHI2)) CHI2 = 1e99;
        cout << "invls: " << x[0] << " g1: " << x[1] << " g2: " << x[2] << " m1: " << x[3] << endl;
        cout << "chi2: " << CHI2 << endl;
        // Return chi2
        return CHI2;
    };

    vector<double> xopt = optimFunction(x_guess, func, 10);
    // Compute best chi2
    // gluon pars correspond to xopt(0)
    gluon_pars = {xopt[0], 0, 0, 0, 0, 0, 0};
    // meson pars are the same as gluon pars
    meson_pars = gluon_pars;
    vector<vector<double> > kernels_pars = {gluon_pars, meson_pars};
    // GNs correspond to x(1), x(2), x(3)
    GNs = {xopt[1], xopt[2], xopt[3]};
    // Update the spectrum
    hqcdp.computeSpectrum(kernels_pars);
    // Update the GNs
    hqcdp.setGNs(GNs);
    // Compute chi2
    double CHI2 = hqcdp.chi2();
    // Print the results
    cout << "Best chi2 found for:" << endl;
    cout << "invls:\t" << gluon_pars[0] << endl;
    // Print the GNs
    cout << "g1:\t" << GNs[0] << endl;
    cout << "g2:\t" << GNs[1] << endl;
    cout << "m1:\t" << GNs[2] << endl;
    cout << "chi2:\t"  << CHI2 << endl;
    const int ndof = hqcdp.NumberOfDegreesOfFreedom();
    cout << "Number of degrees of freedom (Ndof):\t" << ndof << endl;
    cout << "chi2/Ndof:\t" << CHI2/ndof << endl;

    return 0;
}

