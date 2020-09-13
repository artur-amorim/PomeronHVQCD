#include <iostream>
#include <string>
#include <vector>
#include "HolographicVQCD.h"
#include "HQCDP.h"
#include "GluonKernel.h"
#include "MesonKernel.h"
#include "SigmaGammaGamma.h"
#include "methods/optimization/NelderMead.hpp"
#include "schrodinger/chebspec.h"

using namespace std;

int main(int argc, char ** argv)
{
    string data_path;
    double g1, g2, m1;
    double bg = -10.6221, bm = -5.58333;
    if(argc < 5)
    {
        data_path = "expdata/SigmaGammaGamma/SigmaGammaGamma_PDG_data_W_gt_4.txt";
        // Default values for N = 400
        g1 = 1.9765e-7; g2 = 3.5005e-6; m1 = 3.51578e-6;
    }
    else
    {
        data_path = argv[1]; g1 = stod(argv[2]); g2 = stod(argv[3]); m1 = stod(argv[4]);
    }
    cout << "Starting the fit with" << endl;
    cout << "g1: " << g1 << " g2: " << g2 << " m1: " << m1 << endl;

    double mq = hvqcd().QuarkMass();
    SigmaGammaGamma sigma_gg(data_path);

    // Setup Gluon and Meson Kernels and GNs vector
    GluonKernel gluon(2, {0, 0, bg, 0, 0, 0, 0});
    MesonKernel meson(1, {0, 0, bm, 0, 0, 0, 0});
    vector<double> GNs = {g1, g2, m1};

    // Setup HQCDP object
    HQCDP hqcdp;
    hqcdp.addProcessObservable(sigma_gg);
    hqcdp.addKernel(gluon);
    hqcdp.addKernel(meson);
    hqcdp.setGNs(GNs);
    // Compute the spectrum
    // initialise Chebyschev matrices
    chebSetN(1000);
    hqcdp.computeSpectrum();


    vector<double> x_guess = {g1, g2, m1};

       // Definition of the function we want to fit
    auto func = [&hqcdp] (const vector<double> &x)
    {
        // GNs correspond to x(0), x(1), x(2)
        vector<double> gpars = {x[0], x[1], x[2]};
        // Update the GNs
        hqcdp.setGNs(gpars);
        // Compute chi2
        double CHI2 = hqcdp.chi2();
        if(std::isnan(CHI2)) CHI2 = 1e99;
        cout << "g1: " << x[0] << " g2: " << x[1] << " m1: " << x[2] << endl;
        cout << "chi2: " << CHI2 << endl;
        // Return chi2
        return CHI2;
    };

    vector<double> xopt = optimFunction(x_guess, func, 100, 1e-12);
    // Compute best chi2
    GNs = {xopt[0], xopt[1], xopt[2]};
    // Update the GNs
    hqcdp.setGNs(GNs);
    // Compute chi2
    double CHI2 = hqcdp.chi2();
    // Print the results
    cout << "Best chi2 found for:" << endl;
    cout << "g1:\t" << GNs[0] << endl;
    cout << "g2:\t" << GNs[1] << endl;
    cout << "m1:\t" << GNs[2] << endl;
    cout << "chi2:\t"  << CHI2 << endl;
    const int ndof = sigma_gg.expKinematics()[0].size() - GNs.size();
    cout << "Number of degrees of freedom (Ndof):\t" << ndof << endl;
    cout << "chi2/Ndof:\t" << CHI2/ndof << endl;

    return 0;
}