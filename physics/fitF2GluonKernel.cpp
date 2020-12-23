#include <iostream>
#include <string>
#include <vector>
#include "HolographicVQCD.h"
#include "HQCDP.h"
#include "GluonKernel.h"
#include "F2.h"
#include "methods/optimization/NelderMead.hpp"
#include "schrodinger/chebspec.h"

using namespace std;

int main(int argc, char ** argv)
{
    string f2_data_path;
    double ag, bg, cg, dg, eg;
    double g1, g2, g3, g4;
    int N;
    if(argc < 12)
    {
        f2_data_path = "expdata/DIS/F2_data_Q2max_10.txt";
        // Default values for N = 400
        //invls = 0.0439659; ag = 0.305626; bg = -10.1539; cg = 0.822607; dg = -4.85037; eg = -0.251665; fg = 10.0599;
        g1 = -0.213955; g2 = 1.55028; g3 = -3.24117; g4 = 15.9386;
        g1 = 0; g2 = 0; g3 = 0; g4 = 0;
        N = 400;
    }
    else
    {
        f2_data_path = argv[1];
        ag = stod(argv[2]); bg = stod(argv[3]); cg = stod(argv[4]); dg = stod(argv[5]); eg = stod(argv[6]);
        g1 = stod(argv[7]); g2 = stod(argv[8]); g3 = stod(argv[9]); g4 = stod(argv[10]);
        N = stoi(argv[11]);
    }
    cout << "Starting the fit with" << endl;
    cout << "ag: " << ag << " bg: " << bg << " cg: " << cg << " dg: " << dg << " eg: " << eg << endl;
    cout << "g1: " << g1 << " g2: " << g2 << " g3: " << g3 << " g4: " << g4 << endl;

    double mq = hvqcd().QuarkMass();
    F2 f2(f2_data_path);

    // Setup Gluon and Meson Kernels and GNs vector
    GluonKernel gluon(4, {0, ag, bg, cg, dg, eg, 0});
    vector<double> GNs = {g1, g2, g3, g4};

    // Setup HQCDP object
    HQCDP hqcdp;
    hqcdp.addProcessObservable(f2);
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
        cout << "ag = " << x[0] << " bg = " << x[1] << " cg = " << x[2] << " dg = " << x[3] << " eg = " << x[4] << endl;
        cout << "g1 = " << x[5] << " g2 = " << x[6] << " g3 = " << x[7] << " g4 = " << x[8] << endl;
        cout << "chi2: " << CHI2 << endl;
        //double j_sp = hqcdp.getSpectrum()[0].getReggeons()[1].getJ();
        //cout << "Soft Pomeron intercept: " << j_sp << endl;
        //CHI2 +=  ( fabs(j_sp - 1.08) > 0.01 ) ? exp( 1000 * fabs(j_sp - 1.08)) : 0 ;
        //cout << "chi2 + Soft Pomeron constraint: " << CHI2 << endl;
        return CHI2;
    };
    vector<double> deltas = {10, 10, 10, 10, 10, 10, 10, 10, 10};
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
    const int ndof = f2.expKinematics()[0].size() - xopt.size();
    cout << "Number of degrees of freedom (Ndof):\t" << ndof << endl;
    cout << "chi2/Ndof:\t" << CHI2/ndof << endl;

    return 0;
}