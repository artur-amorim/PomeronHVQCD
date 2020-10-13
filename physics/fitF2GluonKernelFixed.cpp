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
    double invls, ag, bg, cg, dg, eg, fg;
    //invls = 0; ag = 0.774954; bg = -12.7148; cg = 3.52521; dg = -4.37674; eg = -0.997324; fg = 0;
    invls = 0; ag = 0.788065; bg = -12.8544; cg = 3.62336; dg = -4.31717; eg = -1.02091; fg = 0;
    double g1, g2, g3, g4;
    int N;
    if(argc < 7)
    {
        f2_data_path = "expdata/DIS/F2_data_Q2max_10.txt";
        // Default values for N = 400
        g1 = -7.38451; g2 = 19.7027; g3 = -0.152453; g4 = 11.5501;
        N = 400;
    }
    else
    {
        f2_data_path = argv[1];
        g1 = stod(argv[2]); g2 = stod(argv[3]); g3 = stod(argv[4]); g4 = stod(argv[5]);
        N = stoi(argv[6]);
    }
    cout << "Starting the fit with" << endl;
    cout << "invls: " << invls << " ag: " << ag << " bg: " << bg << " cg: " << cg << " dg: " << dg << " eg: " << eg << " fg: " << fg << endl;
    cout << "g1: " << g1 << " g2: " << g2 << " g3: " << g3 << " g4: " << g4 << endl;

    double mq = hvqcd().QuarkMass();
    F2 f2(f2_data_path);

    // Setup Gluon and Meson Kernels and GNs vector
    GluonKernel gluon(4, {invls, ag, bg, cg, dg, eg, fg});
    vector<double> GNs = {g1, g2, g3, g4};

    // Setup HQCDP object
    HQCDP hqcdp;
    hqcdp.addProcessObservable(f2);
    hqcdp.addKernel(gluon);
    hqcdp.setGNs(GNs);
    // Compute the spectrum
    // initialise Chebyschev matrices
    chebSetN(N);

    hqcdp.computeSpectrum();

    vector<Reggeon> reggeons = hqcdp.getSpectrum()[0].getReggeons();
    for(int i = 0; i < reggeons.size(); i++) cout << "j" << i << " = " << reggeons[i].getJ() << endl;

    vector<double> x_guess = {g1, g2, g3, g4};

       // Definition of the function we want to fit
    auto func = [&hqcdp] (const vector<double> &x)
    {
        // GNs correspond to x(0), x(1), x(2), x(3)
        vector<double> gpars = {x[0], x[1], x[2], x[3]};
        // Update the GNs
        hqcdp.setGNs(gpars);
        // Compute chi2
        double CHI2 = hqcdp.chi2();
        if(std::isnan(CHI2)) return 1e99;
        cout << "g1 = " << x[0] << " g2 = " << x[1] << " g3 = " << x[2] << " g4 = " << x[3] << endl;
        cout << "chi2: " << CHI2 << endl;
        return CHI2;
    };
    vector<double> deltas = {10, 10, 10, 10};
    vector<double> xopt = optimFunction(x_guess, func, deltas, 1e-12);
    // Compute best chi2
    
    hqcdp.setGNs({xopt[0], xopt[1], xopt[2], xopt[3]});
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