#include <iostream>
#include <string>
#include <vector>
#include "HolographicVQCD.h"
#include "HQCDP.h"
#include "GluonKernel.h"
#include "FL.h"
#include "SigmaGammaProton.h"
#include "methods/optimization/NelderMead.hpp"
#include "schrodinger/chebspec.h"

using namespace std;

int main(int argc, char ** argv)
{
    string fl_data_path, sigma_gp_data_path;
    double invls, ag, bg, cg, dg, eg, fg;
    double g1, g2, g3, g4;
    int N;
    if(argc < 15)
    {
        fl_data_path = "expdata/DIS/FL_data_Q2max_10.txt";
        sigma_gp_data_path = "expdata/SigmaGammaProton/SigmaGammaP_PDG_data_W_gt_461.txt";
        // Default values for N = 400
        invls = 0, ag = 0.0676887; bg = -10.3607; cg = 0.374505; dg = -1.64394; eg = 2.36996; fg = 0;
        g1 = 0.00669477; g2 = -0.0298694; g3 = 0.0541809; g4 = 0.737053;
        N = 400;
    }
    else
    {
        fl_data_path = argv[1]; sigma_gp_data_path = argv[2];
        invls = stod(argv[3]); ag = stod(argv[4]); bg = stod(argv[5]); cg = stod(argv[6]); dg = stod(argv[7]); eg = stod(argv[8]); fg = stod(argv[9]);
        g1 = stod(argv[10]); g2 = stod(argv[11]); g3 = stod(argv[12]); g4 = stod(argv[13]);
        N = stoi(argv[14]);
    }
    cout << "Starting the fit with" << endl;
    cout << "invls: " << invls << " ag: " << ag << " bg: " << bg << " cg: " << cg << " dg: " << dg << " eg: " << eg << " fg: " << fg << endl;
    cout << "g1: " << g1 << " g2: " << g2 << " g3: " << g3 << " g4: " << g4 << endl;

    double mq = hvqcd().QuarkMass();
    FL fl(fl_data_path);
    SigmaGammaProton sigma_gp(sigma_gp_data_path);

    // Setup Gluon and Meson Kernels and GNs vector
    GluonKernel gluon(4, {invls, ag, bg, cg, dg, eg, fg});
    vector<double> GNs = {g1, g2, g3, g4};

    // Setup HQCDP object
    HQCDP hqcdp;
    hqcdp.addProcessObservable(fl);
    hqcdp.addProcessObservable(sigma_gp);
    hqcdp.addKernel(gluon);
    hqcdp.setGNs(GNs);
    // initialise Chebyschev matrices
    chebSetN(N);

    vector<double> x_guess = {ag, bg, cg, dg, eg, g1, g2, g3, g4};

    // Definition of the function we want to fit
    auto func = [&hqcdp] (const vector<double> &x)
    {
        // Compute the spectrum with the gluon kernel parameters
        hqcdp.computeSpectrum({{0, x[0], x[1], x[2], x[3], x[4], 0}});
        // GNs correspond to x(7), x(8), x(9), x(10)
        vector<double> gpars = {x[5], x[6], x[7], x[8]};
        // Update the GNs
        hqcdp.setGNs(gpars);
        // Compute chi2
        double CHI2 = hqcdp.chi2();
        if(std::isnan(CHI2)) return 1e99;
        cout << "invls: " << 0 << " ag: " << x[0] << " bg: " << x[1] << " cg: " << x[2] << " dg: " << x[3] << " eg: " << x[4] << " fg: " << 0 << endl;
        cout << "g1: " << x[5] << " g2: " << x[6] << " g3: " << x[7] << " g4: " << x[8] << endl;
        cout << "chi2: " << CHI2 << endl;
        return CHI2;
    };
    vector<double> deltas = {10, 10, 10, 10, 10, 1, 1, 1, 1};
    vector<double> xopt = optimFunction(x_guess, func, deltas, 1e-12);
    
    // Compute best chi2
    hqcdp.computeSpectrum({{0, xopt[0], xopt[1], xopt[2], xopt[3], xopt[4], 0}});
    hqcdp.setGNs({xopt[5], xopt[6], xopt[7], xopt[8]});
    // Compute chi2
    double CHI2 = hqcdp.chi2();
    // Print the results
    cout << "Best chi2 found for: " << endl;
    for(int i = 0; i < xopt.size(); i++) cout << xopt[i] << '\t';
    cout << endl;
    cout << "chi2:\t"  << CHI2 << endl;
    const int ndof = fl.expKinematics()[0].size() + sigma_gp.expKinematics()[0].size() - xopt.size();
    cout << "Number of degrees of freedom (Ndof):\t" << ndof << endl;
    cout << "chi2/Ndof:\t" << CHI2/ndof << endl;

    return 0;
}