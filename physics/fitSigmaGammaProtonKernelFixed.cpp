#include <iostream>
#include <string>
#include <vector>
#include "HolographicVQCD.h"
#include "HQCDP.h"
#include "GluonKernel.h"
#include "MesonKernel.h"
#include "SigmaGammaProton.h"
#include "methods/optimization/NelderMead.hpp"
#include "schrodinger/chebspec.h"

using namespace std;

int main(int argc, char ** argv)
{
    string data_path;
    double g1, g2, m1;
    //double invls_g = 0.253093, ag = -0.322648, bg = -0.0559966, cg = -0.7771, dg = 0.0145733, eg = 0.256595, fg = 0.131464;
    //double invls_m = -0.00513779, am = 3.57895, bm = 0.0192416, cm = 0.486709, dm = -0.0102358, em = 6.43404, fm = 0.333183;
    double bg = -10.6221, bm = -5.58333;
    if(argc < 5)
    {
        data_path = "expdata/SigmaGammaProton/SigmaGammaP_PDG_data_W_gt_461.txt";
        g1 = 0; g2 = 0; m1 = 0;
    }
    else
    {
        data_path = argv[1]; g1 = stod(argv[2]); g2 = stod(argv[3]); m1 = stod(argv[4]);
    }
    cout << "Starting the fit with" << endl;
    cout << "g1: " << g1 << " g2: " << g2 << " m1: " << m1 << endl;

    double mq = hvqcd().QuarkMass();
    SigmaGammaProton sigma_gp(data_path);

    // Setup Gluon and Meson Kernels and GNs vector
    GluonKernel gluon(2, {0, 0, bg, 0, 0, 0, 0});
    MesonKernel meson(1, {0, 0, bm, 0, 0, 0, 0});
    vector<double> GNs = {g1, g2, m1};

    // Setup HQCDP object
    HQCDP hqcdp;
    hqcdp.addProcessObservable(sigma_gp);
    hqcdp.addKernel(gluon);
    hqcdp.addKernel(meson);
    hqcdp.setGNs(GNs);
    // Compute the spectrum
    // initialise Chebyschev matrices
    chebSetN(400);
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
        else CHI2 += exp(-x[0]) + exp(-x[1]) + exp(-x[2]);
        cout << "g1: " << x[0] << " g2: " << x[1] << " m1: " << x[2] << endl;
        cout << "chi2: " << CHI2 << endl;
        // Return chi2
        return CHI2;
    };

    vector<double> xopt = optimFunction(x_guess, func, 10, 1e-12);
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
    const int ndof = sigma_gp.expKinematics()[0].size() - GNs.size();
    cout << "Number of degrees of freedom (Ndof):\t" << ndof << endl;
    cout << "chi2/Ndof:\t" << CHI2/ndof << endl;

    return 0;
}