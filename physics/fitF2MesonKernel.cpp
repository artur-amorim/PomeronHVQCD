#include <iostream>
#include <vector>
#include "HolographicVQCD.h"
#include "F2.h"
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
    double g1g, g2g, g3g, g4g, g1m;
    if (argc < 7)
    {
        invls = 0; ag = -0.788065; bg = 12.8544; cg = -3.62336; dg = 4.31717; eg = 1.02091; fg = 0;
        am = 0; bm = 5.58333; cm = 0; dm = 0; em = 0; fm = 0;
        g1g = -7.38451; g2g = 19.7027; g3g = -0.152453; g4g = 11.5501; g1m = 1;
    }
    else
    {
        am = stod(argv[1]); bm = stod(argv[2]); cm = stod(argv[3]); dm = stod(argv[4]);
        em = stod(argv[5]); g1m = stod(argv[6]);
    }

    cout << "Starting the fit with" << endl;
    cout << "invls: " << invls << " ag: " << ag << " bg: " << bg << " cg: " << cg << " dg: " << dg;
    cout << " eg: " << eg <<  " fg: " << fg << endl;
    cout << "am: " << am << " bm: " << bm << " cm: " << cm << " dm: " << dm << " em: " << em <<  " fm: " << fm << endl;
    cout << "g1g: " << g1g << " g2g: " << g2g << " g3g: " << g3g << " g4g: " << g4g << endl;
    cout << "g1m: " << g1m << endl;

    double mq = hvqcd().QuarkMass();
    F2 f2("expdata/DIS/F2_data_Q2max_10.txt");

    // Setup Gluon and Meson Kernels and GNs vector
    vector<double> gluon_pars = {invls, ag, bg, cg, dg, eg, fg};
    GluonKernel gluon(4, gluon_pars);
    vector<double> meson_pars = {invls, am, bm, cm, dm, em, fm};
    MesonKernel meson(1, meson_pars);
    vector<double> GNs = {g1g, g2g, g3g, g4g, g1m};

    // Setup HQCDP object
    HQCDP hqcdp;
    hqcdp.addProcessObservable(f2);
    hqcdp.addKernel(gluon);
    hqcdp.addKernel(meson);
    hqcdp.setGNs(GNs);
    
    // initialise Chebyschev matrices
    chebSetN(400);

    vector<double> x_guess = {am, bm, cm, dm, em, g1m};

    // Definition of the function we want to fit
    auto func = [&hqcdp, &gluon_pars, &GNs] (const vector<double> &x)
    {
        // meson pars correspond to x(0), x(7), x(8), x(9), x(10), x(11) and x(12)
        vector<double> meson_pars = {0, x[0], x[1], x[2], x[3], x[4], 0};
        vector<vector<double> > kernels_pars = {gluon_pars, meson_pars};
        // GNs correspond to x(13), x(14), x(15)
        vector<double> gpars = {GNs[0], GNs[1], GNs[2], GNs[3], x[5]};
        // Update the spectrum
        hqcdp.computeSpectrum(kernels_pars);
        // Update the GNs
        hqcdp.setGNs(gpars);
        // Compute chi2
        double CHI2 = hqcdp.chi2();
        if(std::isnan(CHI2)) CHI2 = 1e99;
        cout << "am = " << x[0] << " bm = " << x[1] << " cm = " << x[2] << " dm = " << x[3] << " em = " << x[4] <<  " fm = " << 0 << endl;
        cout << "g1m = " << x[5] << endl;
        cout << "chi2: " << CHI2 << endl;
        // Return chi2
        return CHI2;
    };

    vector<double> deltas = {10, 10, 10, 10, 10, 1};
    vector<double> xopt = optimFunction(x_guess, func, deltas);

    // Compute best chi2
    meson_pars = {0, xopt[0], xopt[1], xopt[2], xopt[3], xopt[4], 0};
    vector<vector<double> > kernels_pars = {gluon_pars, meson_pars};
    // GNs correspond to x(13), x(14), x(15)
    GNs = {GNs[0], GNs[1], GNs[2], GNs[3], xopt[5]};
    // Update the spectrum
    hqcdp.computeSpectrum(kernels_pars);
    // Update the GNs
    hqcdp.setGNs(GNs);
    // Compute chi2
    double CHI2 = hqcdp.chi2();
    // Print the results
    cout << "Best chi2 found for:" << endl;
    cout << "invls:\t" << gluon_pars[0] << endl;
    cout << "ag:\t"     << gluon_pars[1] << endl;
    cout << "bg:\t"     << gluon_pars[2] << endl;
    cout << "cg:\t"     << gluon_pars[3] << endl;
    cout << "dg:\t"     << gluon_pars[4] << endl;
    cout << "eg:\t"     << gluon_pars[5] << endl;
    cout << "fg:\t"     << gluon_pars[6] << endl;
    cout << "am:\t"     << meson_pars[1] << endl;
    cout << "bm:\t"     << meson_pars[2] << endl;
    cout << "cm:\t"     << meson_pars[3] << endl;
    cout << "dm:\t"     << meson_pars[4] << endl;
    cout << "em:\t"     << meson_pars[5] << endl;
    cout << "fm:\t"     << meson_pars[6] << endl;
    // Print the GNs
    cout << "g1g:\t" << GNs[0] << endl;
    cout << "g2g:\t" << GNs[1] << endl;
    cout << "g3g:\t" << GNs[2] << endl;
    cout << "g4g:\t" << GNs[3] << endl;
    cout << "g1m:\t" << GNs[4] << endl;
    cout << "chi2:\t"  << CHI2 << endl;
    const int ndof = hqcdp.NumberOfDegreesOfFreedom();
    cout << "Number of degrees of freedom (Ndof):\t" << ndof << endl;
    cout << "chi2/Ndof:\t" << CHI2/ndof << endl;

    return 0;
}