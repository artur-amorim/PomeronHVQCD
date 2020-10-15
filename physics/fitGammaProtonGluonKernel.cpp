#include <iostream>
#include <string>
#include <vector>
#include "HolographicVQCD.h"
#include "HQCDP.h"
#include "GluonKernel.h"
#include "F2.h"
#include "FL.h"
#include "SigmaGammaProton.h"
#include "methods/optimization/NelderMead.hpp"
#include "schrodinger/chebspec.h"

using namespace std;

int main(int argc, char ** argv)
{
    string f2_data_path, fl_data_path, sigma_gp_data_path;
    double invls, ag, bg, cg, dg, eg, fg;
    double g1, g2, g3, g4;
    int N;
    if(argc < 16)
    {
        f2_data_path = "expdata/DIS/F2_data_Q2max_10.txt";
        fl_data_path = "expdata/DIS/FL_data_Q2max_10.txt";
        sigma_gp_data_path = "expdata/SigmaGammaProton/SigmaGammaP_PDG_data_W_gt_461.txt";
        // Default values for N = 400
        //invls = 0.0439659; ag = 0.305626; bg = -10.1539; cg = 0.822607; dg = -4.85037; eg = -0.251665; fg = 10.0599;
        invls = 0; ag = 0.0486791; bg = -11.8035; cg = 1.28621; dg = -4.8624; eg = -0.366335; fg = 0;
        g1 = 0.0220324; g2 = -0.0174865; g3 = 0.148405; g4 = -38.6718;
        N = 400;
    }
    else
    {
        f2_data_path = argv[1]; fl_data_path = argv[2]; sigma_gp_data_path = argv[3];
        invls = stod(argv[4]); ag = stod(argv[5]); bg = stod(argv[6]); cg = stod(argv[7]); dg = stod(argv[8]); eg = stod(argv[9]); fg = stod(argv[10]);
        g1 = stod(argv[11]); g2 = stod(argv[12]); g3 = stod(argv[13]); g4 = stod(argv[14]);
        N = stoi(argv[15]);
    }
    cout << "Starting the fit with" << endl;
    cout << "invls: " << invls << " ag: " << ag << " bg: " << bg << " cg: " << cg << " dg: " << dg << " eg: " << eg << " fg: " << fg << endl;
    cout << "g1: " << g1 << " g2: " << g2 << " g3: " << g3 << " g4: " << g4 << endl;

    double mq = hvqcd().QuarkMass();
    F2 f2(f2_data_path);
    FL fl(fl_data_path);
    SigmaGammaProton sigma_gp(sigma_gp_data_path);

    // Setup Gluon and Meson Kernels and GNs vector
    GluonKernel gluon(4, {invls, ag, bg, cg, dg, eg, fg});
    vector<double> GNs = {g1, g2, g3, g4};

    // Setup HQCDP object
    HQCDP hqcdp;
    hqcdp.addProcessObservable(f2);
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
        cout << "ag: " << x[0] << " bg: " << x[1] << " cg: " << x[2] << " dg: " << x[3] << " eg: " << x[4] << endl;
        cout << "g1: " << x[5] << " g2: " << x[6] << " g3: " << x[7] << " g4: " << x[8] << endl;
        cout << "chi2: " << CHI2 << endl;
        double j_sp = hqcdp.getSpectrum()[0].getReggeons()[1].getJ();
        cout << "Soft Pomeron intercept: " << j_sp << endl;
        CHI2 += 10000 * std::fabs(j_sp - 1.08);
        cout << "chi2 + Soft Pomeron constraint: " << CHI2 << endl;
        return CHI2;
    };
    vector<double> deltas = {-1, -1, -5, -1, -1, 1, 1, 1, 10};
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
    const int ndof = f2.expKinematics()[0].size() + fl.expKinematics()[0].size() + sigma_gp.expKinematics()[0].size() - xopt.size();
    cout << "Number of degrees of freedom (Ndof):\t" << ndof << endl;
    cout << "chi2/Ndof:\t" << CHI2/ndof << endl;

    return 0;
}