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
    double g1g, g2g, g1m;
    if (argc < 17)
    {
        invls = 0.254119; ag = 13.9538; bg = 0.921665; cg = 2.03904; dg = -2.7305; eg = -0.473787;
        fg = -0.517072;
        am = 0; bm = 0; cm = 0, dm = 0; em = 0; fm = 0;
        g1g = 1; g2g = 1; g1m = 1;
    }
    else
    {
        invls = stod(argv[1]); ag = stod(argv[2]); bg = stod(argv[3]); cg = stod(argv[4]); dg = stod(argv[5]);
        eg = stod(argv[6]); fg = stod(argv[7]);
        am = stod(argv[8]); bm = stod(argv[9]); cm = stod(argv[10]); dm = stod(argv[11]);
        em = stod(argv[12]); fm = stod(argv[13]);
        g1g = stod(argv[14]); g2g = stod(argv[15]); g1m = stod(argv[16]);
    }

    cout << "Starting the fit with" << endl;
    cout << "invls: " << invls << " ag: " << ag << " bg: " << bg << " cg: " << cg << " dg: " << dg;
    cout << " eg: " << eg <<  " fg: " << fg << endl;
    cout << "am: " << am << " bm: " << bm << " cm: " << cm << " dm: " << dm << " em: " << em <<  " fm: " << fm << endl;
    cout << "g1g: " << g1g << " g2g: " << g2g << " g1m: " << g1m << endl;

    double mq = hvqcd().QuarkMass();
    F2 f2("expdata/F2_data.txt");
    FL fl("expdata/FL_data.txt");

    // Setup Gluon and Meson Kernels and GNs vector
    vector<double> gluon_pars = {invls, ag, bg, cg, dg, eg, fg};
    GluonKernel gluon(2, gluon_pars);
    vector<double> meson_pars = {invls, am, bm, cm, dm, em, fm};
    MesonKernel meson(1, meson_pars);
    vector<double> GNs = {g1g, g2g, g1m};

    // Setup HQCDP object
    HQCDP hqcdp(false ,0);
    hqcdp.addProcessObservable(f2);
    hqcdp.addProcessObservable(fl);
    hqcdp.addKernel(gluon);
    hqcdp.addKernel(meson);
    hqcdp.setGNs(GNs);
    
    // initialise Chebyschev matrices
    chebSetN(1000);

    vector<double> x_guess = {invls, ag, bg, cg, dg, eg, fg, am, bm, cm, dm, em, fm, g1g, g2g, g1m};

    // Definition of the function we want to fit
    auto func = [&hqcdp] (const vector<double> &x)
    {
        // gluon pars correspond to x(0), x(1), x(2), x(3), x(4), x(5) and x(6)
        vector<double> gluon_pars = {x[0], x[1], x[2], x[3], x[4], x[5], x[6]};
        // meson pars correspond to x(0), x(7), x(8), x(9), x(10), x(11) and x(12)
        vector<double> meson_pars = {x[0], x[7], x[8], x[9], x[10], x[11], x[12]};
        vector<vector<double> > kernels_pars = {gluon_pars, meson_pars};
        // GNs correspond to x(13), x(14), x(15)
        vector<double> gpars = {x[13], x[14], x[15]};
        // Update the spectrum
        hqcdp.computeSpectrum(kernels_pars);
        // Update the GNs
        hqcdp.setGNs(gpars);
        // Compute chi2
        double CHI2 = hqcdp.chi2();
        if(std::isnan(CHI2)) CHI2 = 1e99;
        cout << "invls: " << x[0] << " ag: " << x[1] << " bg: " << x[2] << " cg: " << x[3] << " dg: " << x[4];
        cout << " eg: " << x[5] <<  " fg: " << x[6] << endl;
        cout << "am: " << x[7] << " bm: " << x[8] << " cm: " << x[9] << " dm: " << x[10] << " em: " << x[11] <<  " fm: " << x[12] << endl;
        cout << "g1g: " << x[13] << " g2g: " << x[14] << " g1m: " << x[15] << endl;
        cout << "chi2: " << CHI2 << endl;
        // Return chi2
        return CHI2;
    };

    vector<double> xopt = optimFunction(x_guess, func, 10);

    // Compute best chi2
    // gluon pars correspond to xopt(0), xopt(1), xopt(2), xopt(3), xopt(4), xopt(5), xopt(6)
    gluon_pars = {xopt[0], xopt[1], xopt[2], xopt[3], xopt[4], xopt[5], xopt[6]};
    // meson pars correspond to xopt(0), xopt(7), xopt(8), xopt(9), xopt(10), xopt(11) and xopt(12)
    meson_pars = {xopt[0], xopt[7], xopt[8], xopt[9], xopt[10], xopt[11], xopt[12]};
    vector<vector<double> > kernels_pars = {gluon_pars, meson_pars};
    // GNs correspond to x(13), x(14), x(15)
    GNs = {xopt[13], xopt[14], xopt[15]};
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
    cout << "g1m:\t" << GNs[2] << endl;
    cout << "chi2:\t"  << CHI2 << endl;
    const int ndof = hqcdp.NumberOfDegreesOfFreedom();
    cout << "Number of degrees of freedom (Ndof):\t" << ndof << endl;
    cout << "chi2/Ndof:\t" << CHI2/ndof << endl;

    return 0;
}