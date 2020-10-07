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
    double ag, am;
    double bg, bm;
    double cg, cm;
    double dg, dm;
    double eg, em;
    double g1, g2, m1, m2;
    int N;
    if(argc < 18)
    {
        fl_data_path = "expdata/DIS/FL_data_Q2max_10.txt";
        sigma_gp_data_path = "expdata/SigmaGammaProton/SigmaGammaP_PDG_data_W_gt_461.txt";
        // Default values for N = 400
        ag = 0, bg = -10.6221, cg = 0, dg = 0, eg = 0;
        am = 0, bm = -5.58333, cm = 0, dm = 0, em = 0;
        g1 = 0.00907812; g2 = 0.0383623; m1 = 0.0272387, m2 = 0;
        N = 400;
    }
    else
    {
        fl_data_path = argv[1]; sigma_gp_data_path = argv[2];
        ag = stod(argv[3]); bg = stod(argv[4]); cg = stod(argv[5]); dg = stod(argv[6]); eg = stod(argv[7]);
        am = stod(argv[8]); bm = stod(argv[9]); cm = stod(argv[10]); dm = stod(argv[11]); em = stod(argv[12]);
        g1 = stod(argv[13]); g2 = stod(argv[14]); m1 = stod(argv[15]); m2 = stod(argv[16]);
        N = stoi(argv[17]);
    }
    cout << "Starting the fit with" << endl;
    cout << "ag: " << ag << " bg: " << bg << " cg: " << cg << " dg: " << dg << " eg: " << eg << endl;
    cout << "am: " << am << " bm: " << bm << " cm: " << cm << " dm: " << dm << " em: " << em << endl;
    cout << "g1: " << g1 << " g2: " << g2 << " m1: " << m1 << " m2: " << m2 << endl;

    double mq = hvqcd().QuarkMass();
    FL fl(fl_data_path);
    SigmaGammaProton sigma_gp(sigma_gp_data_path);

    // Setup Gluon and Meson Kernels and GNs vector
    GluonKernel gluon(2, {0, ag, bg, cg, dg, eg, 0});
    MesonKernel meson(2, {0, am, bm, cm, dm, em, 0});
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
    chebSetN(N);

    vector<double> x_guess = {ag, bg, cg, dg, eg, am, bm, cm, dm, em, g1, g2, m1, m2};

       // Definition of the function we want to fit
    auto func = [&hqcdp] (const vector<double> &x)
    {
        // gluon pars correspond to x(0)
        vector<double> gluon_pars = {0, x[0], x[1], x[2], x[3], x[4], 0};
        // meson pars correspond are the same as gluon_pars
        vector<double> meson_pars = {0, x[5], x[6], x[7], x[8], x[9], 0};;
        vector<vector<double> > kernels_pars = {gluon_pars, meson_pars};
        // GNs correspond to x(1), x(2), x(3)
        vector<double> gpars = {x[10], x[11], x[12], x[13]};
        // Update the spectrum
        hqcdp.computeSpectrum(kernels_pars);
        // Update the GNs
        hqcdp.setGNs(gpars);
        // Compute chi2
        double CHI2 = hqcdp.chi2();
        if(std::isnan(CHI2)) return 1e99;
        cout << "chi2: " << CHI2 << endl;
        /*vector<Reggeon> reggeons = hqcdp.getSpectrum()[0].getReggeons();
        double soft_pomeron_intercept = reggeons[1].getJ();
        double meson_intercept = reggeons[2].getJ();
        // Add soft pomeron and meson constraints
        CHI2 += 10000 * fabs(soft_pomeron_intercept - 1.0808);
        CHI2 += 10000 * fabs(meson_intercept - 0.5475);*/
        cout << "ag: " << x[0] << " bg: " << x[1] << " cg: " << x[2] << " dg: " << x[3] << " eg: " << x[4] << endl;
        cout << "am: " << x[5] << " bm: " << x[6] << " cm: " << x[7] << " dm: " << x[8] << " em: " << x[9] << endl;
        cout << "g1: " << x[10] << " g2: " << x[11] << " m1: " << x[12] << " m2: " << x[13] << endl;
        /*cout << "Soft Pomeron Intercept: " << soft_pomeron_intercept << endl;
        cout << "Meson Intercept: " << meson_intercept << endl;*/
        cout << "chi2: " << CHI2 << endl;
        // Return chi2
        return CHI2;
    };
    vector<double> deltas = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 1, 1, 1, 1};
    vector<double> xopt = optimFunction(x_guess, func, deltas, 1e-12);
    // Compute best chi2
    // Update the GNs
    hqcdp.computeSpectrum({{0, xopt[0], xopt[1], xopt[2], xopt[3], xopt[4], 0}, {0, xopt[5], xopt[6], xopt[7], xopt[8], xopt[9], 0}});
    hqcdp.setGNs({xopt[10], xopt[11], xopt[12], xopt[13]});
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