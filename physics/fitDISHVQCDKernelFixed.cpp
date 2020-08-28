#include <iostream>
#include <vector>
#include "HolographicVQCD.h"
#include "F2.h"
#include "FL.h"
#include "GluonKernel.h"
#include "MesonKernel.h"
#include "HQCDP.h"
#include "schrodinger/chebspec.h"
#include "schrodinger/schrodinger.h"
#include "methods/optimization/NelderMead.hpp"


using namespace std;

int main(int argc, char ** argv)
{
    double invls, ag, bg, cg, dg, eg, fg;
    double am, bm, cm, dm, em, fm;
    invls = 0; ag = 0; bg = -10.6221; cg = 0; dg = 0; eg = 0; fg = 0;
    am = 0; bm = -5.58333; cm = 0, dm = 0; em = 0; fm = 0;
    double g1g, g2g, g1m;
    if (argc < 4)
    {
        g1g = 8.53452e-6; g2g = 0.00206828; g1m = 0.00052165;
    }
    else
    {
        g1g = stod(argv[1]); g2g = stod(argv[2]); g1m = stod(argv[3]);
    }

    cout << "Starting the fit with" << endl;
    cout << "g1g: " << g1g << " g2g: " << g2g << " g1m: " << g1m << endl;

    double mq = hvqcd().QuarkMass();
    F2 f2("expdata/DIS/F2_data.txt");
    FL fl("expdata/DIS/FL_data.txt");

    // Setup Gluon and Meson Kernels and GNs vector
    GluonKernel gluon(2, {0, 0, bg, 0, 0, 0, 0});
    MesonKernel meson(1, {0, 0, bm, 0, 0, 0, 0});
    vector<double> GNs = {g1g, g2g, g1m};

    // Setup HQCDP object
    HQCDP hqcdp;
    hqcdp.addProcessObservable(f2);
    hqcdp.addProcessObservable(fl);
    hqcdp.addKernel(gluon);
    hqcdp.addKernel(meson);
    hqcdp.setGNs(GNs);
    
    // initialise Chebyschev matrices and compute the spectrum
    chebSetN(400);
    // Update the spectrum
    hqcdp.computeSpectrum();

    vector<double> x_guess = {g1g, g2g, g1m};

    // Definition of the function we want to fit
    auto func = [&hqcdp] (const vector<double> &x)
    {
        vector<double> gpars = {x[0], x[1], x[2]};
        // Update the GNs
        hqcdp.setGNs(gpars);
        // Compute chi2
        double CHI2 = hqcdp.chi2();
        if(std::isnan(CHI2)) CHI2 = 1e99;
        cout << "g1g: " << x[0] << " g2g: " << x[1] << " g1m: " << x[2] << endl;
        cout << "chi2: " << CHI2 << endl;
        // Return chi2
        return CHI2;
    };

    vector<double> xopt = optimFunction(x_guess, func, 10, 1e-12);

    // Compute best chi2
    // GNs correspond to x(0), x(1), x(2)
    GNs = {xopt[0], xopt[1], xopt[2]};
    // Update the GNs
    hqcdp.setGNs(GNs);
    // Compute chi2
    double CHI2 = hqcdp.chi2();
    // Print the results
    cout << "Best chi2 found for:" << endl;
    // Print the GNs
    cout << "g1g:\t" << GNs[0] << endl;
    cout << "g2g:\t" << GNs[1] << endl;
    cout << "g1m:\t" << GNs[2] << endl;
    cout << "chi2:\t"  << CHI2 << endl;
    const int ndof = f2.expVal().size() + fl.expVal().size() - GNs.size();
    cout << "Number of degrees of freedom (Ndof):\t" << ndof << endl;
    cout << "chi2/Ndof:\t" << CHI2/ndof << endl;

    return 0;
}