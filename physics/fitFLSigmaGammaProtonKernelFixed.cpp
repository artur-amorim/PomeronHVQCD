#include <iostream>
#include <string>
#include <vector>
#include "HolographicVQCD.h"
#include "HQCDP.h"
#include "GluonKernel.h"
#include "MesonKernel.h"
#include "FL.h"
#include "SigmaGammaProton.h"
#include "methods/optimization/NelderMead.hpp"
#include "schrodinger/chebspec.h"

using namespace std;

int main(int argc, char ** argv)
{
    string fl_data_path, sigma_gp_data_path;
    double invls = 0, fg = 0, fm = 0;
    double ag = 0, bg = -10.6221, cg = 0, dg = 0, eg = 0;
    double am = 0, bm = -5.58333, cm = 0, dm = 0, em = 0;
    double g1, g2, m1, m2;
    if(argc < 7)
    {
        fl_data_path = "expdata/DIS/FL_data_Q2max_10.txt";
        sigma_gp_data_path = "expdata/SigmaGammaProton/SigmaGammaP_PDG_data_W_gt_461.txt";
        // Default values for N = 1000
        g1 = 0.00907812; g2 = 0.0383623; m1 = 0.0272387, m2 = 0;
    }
    else
    {
        fl_data_path = argv[1]; sigma_gp_data_path = argv[2];
        g1 = stod(argv[3]); g2 = stod(argv[4]); m1 = stod(argv[5]); m2 = stod(argv[6]);
    }
    cout << "Starting the fit with" << endl;
    cout << "g1: " << g1 << " g2: " << g2 << " m1: " << m1 << " m2: " << m2 << endl;

    double mq = hvqcd().QuarkMass();
    FL fl(fl_data_path);
    SigmaGammaProton sigma_gp(sigma_gp_data_path);

    // Setup Gluon and Meson Kernels and GNs vector
    GluonKernel gluon(2, {invls, ag, bg, cg, dg, eg, fg});
    MesonKernel meson(2, {invls, am, bm, cm, dm, em, fm});
    vector<double> GNs = {g1, g2, m1, m2};

    // Setup HQCDP object
    HQCDP hqcdp;
    hqcdp.addProcessObservable(fl);
    hqcdp.addProcessObservable(sigma_gp);
    hqcdp.addKernel(gluon);
    hqcdp.addKernel(meson);
    hqcdp.setGNs(GNs);
    // Compute the spectrum
    // initialise Chebyschev matrices
    chebSetN(400);
    hqcdp.computeSpectrum();

    vector<double> x_guess = {g1, g2, m1, m2};

       // Definition of the function we want to fit
    auto func = [&hqcdp] (const vector<double> &x)
    {
        // GNs correspond to x(0), x(1), x(2)
        vector<double> gpars = {x[0], x[1], x[2], x[3]};
        // Update the GNs
        hqcdp.setGNs(gpars);
        // Compute chi2
        double CHI2 = hqcdp.chi2();
        if(std::isnan(CHI2)) CHI2 = 1e99;
        //else CHI2 += exp(-x[0]) + exp(-x[1]) + exp(-x[2]);
        cout << "g1: " << x[0] << " g2: " << x[1] << " m1: " << x[2] << " m2: " << x[3] << endl;
        cout << "chi2: " << CHI2 << endl;
        // Return chi2
        return CHI2;
    };

    vector<double> xopt = optimFunction(x_guess, func, 10, 1e-12);
    // Compute best chi2
    GNs = {xopt[0], xopt[1], xopt[2], xopt[3]};
    // Update the GNs
    hqcdp.setGNs(GNs);
    // Compute chi2
    double CHI2 = hqcdp.chi2();
    // Print the results
    cout << "Best chi2 found for:" << endl;
    cout << "g1:\t" << GNs[0] << endl;
    cout << "g2:\t" << GNs[1] << endl;
    cout << "m1:\t" << GNs[2] << endl;
    cout << "m2:\t" << GNs[3] << endl;
    cout << "chi2:\t"  << CHI2 << endl;
    const int ndof = fl.expKinematics()[0].size() + sigma_gp.expKinematics()[0].size() - GNs.size();
    cout << "Number of degrees of freedom (Ndof):\t" << ndof << endl;
    cout << "chi2/Ndof:\t" << CHI2/ndof << endl;

    return 0;
}