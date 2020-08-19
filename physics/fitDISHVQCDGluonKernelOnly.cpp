#include <iostream>
#include <vector>
#include "HolographicVQCD.h"
#include "F2.h"
#include "FL.h"
#include "GluonKernel.h"
#include "HQCDP.h"
#include "schrodinger/chebspec.h"
#include "schrodinger/schrodinger.h"
#include "methods/optimization/NelderMead.hpp"


using namespace std;

int main(int argc, char ** argv)
{
    double invls, a, b, c, d, e, f;
    double g1, g2, g3, g4;
    if (argc < 12)
    {
        invls = 0.254119; a = 13.9538; b = 0.921665; c = 2.03904; d = -2.7305; e = -0.473787;
        f = -0.517072;
        g1 = 1; g2 = -1; g3 = 1; g4 = -1;
    }
    else
    {
        invls = stod(argv[1]); a = stod(argv[2]); b = stod(argv[3]); c = stod(argv[4]); d= stod(argv[5]);
        e = stod(argv[6]); f = stod(argv[7]);
        g1 = stod(argv[8]); g2 = stod(argv[9]); g3 = stod(argv[10]); g4 = stod(argv[11]);
    }

    cout << "Starting the fit with" << endl;
    cout << "invls: " << invls << " a: " << a << " b: " << b << " c: " << c << " d: " << d;
    cout << " e: " << e <<  " f: " << f << endl;
    cout << " g1: " << g1 << " g2: " << g2 << " g3: " << g3 << " g4: " << g4 << endl;

    double mq = hvqcd().QuarkMass();
    F2 f2("expdata/DIS/F2_data.txt");
    FL fl("expdata/DIS/FL_data.txt");

    // Setup Gluon Kernel and GNs vector
    vector<double> gluon_pars = {invls, a, b, c, d, e, f};
    GluonKernel gluon(4, gluon_pars);
    vector<double> GNs = {g1, g2, g3, g4};

    // Setup HQCDP object
    HQCDP hqcdp(false ,0);
    hqcdp.addProcessObservable(f2);
    hqcdp.addProcessObservable(fl);
    hqcdp.addKernel(gluon);
    hqcdp.setGNs(GNs);
    
    // initialise Chebyschev matrices
    chebSetN(1000);

    vector<double> x_guess = {invls, a, b, c, d, e, f, g1, g2, g3, g4};

    // Definition of the function we want to fit
    auto func = [&hqcdp] (const vector<double> &x)
    {
        // Kernel pars correspond to x(0), x(1), x(2), x(3), x(4), x(5) and x(6)
        vector<double> gluon_pars = {x[0], x[1], x[2], x[3], x[4], x[5], x[6]};
        vector<vector<double> > kernels_pars = {gluon_pars};
        // GNs correspond to x(7), x(8), x(9), x(10)
        vector<double> gpars = {x[7], x[8], x[9], x[10]};
        // Update the spectrum
        hqcdp.computeSpectrum(kernels_pars);
        // Update the GNs
        hqcdp.setGNs(gpars);
        // Compute chi2
        double CHI2 = hqcdp.chi2();
        if(std::isnan(CHI2)) CHI2 = 1e99;
        cout << "invls: " << x[0] << " a: " << x[1] << " b: " << x[2] << " c: " << x[3] << " d: " << x[4];
        cout << " e: " << x[5] << " f: " << x[6];
        cout << " g1: " << x[7] << " g2: " << x[8] << " g3: " << x[9] <<  " g4: " << x[10] << endl;
        cout << "chi2: " << CHI2 << endl;
        // Return chi2
        return CHI2;
    };

    vector<double> xopt = optimFunction(x_guess, func, 10);

    // Compute best chi2
    // Kernel pars correspond to xopt(0), xopt(1), xopt(2), xopt(3), xopt(4), xopt(5), xopt(6)
    gluon_pars = {xopt[0], xopt[1], xopt[2], xopt[3], xopt[4], xopt[5], xopt[6]};
    vector<vector<double> > kernels_pars = {gluon_pars};
    // GNs correspond to x(7), x(8), x(9), x(10)
    GNs = {xopt[7], xopt[8], xopt[9], xopt[10]};
    // Update the spectrum
    hqcdp.computeSpectrum(kernels_pars);
    // Update the GNs
    hqcdp.setGNs(GNs);
    // Compute chi2
    double CHI2 = hqcdp.chi2();
    // Print the results
    cout << "Best chi2 found for:" << endl;
    cout << "invls:\t" << gluon_pars[0] << endl;
    cout << "a:\t"     << gluon_pars[1] << endl;
    cout << "b:\t"     << gluon_pars[2] << endl;
    cout << "c:\t"     << gluon_pars[3] << endl;
    cout << "d:\t"     << gluon_pars[4] << endl;
    cout << "e:\t"     << gluon_pars[5] << endl;
    cout << "f:\t"     << gluon_pars[6] << endl;
    // Print the GNs
    for(int i = 0; i < GNs.size(); i++) cout << "g" << to_string(i+1) << '\t' << GNs[i] << endl;
    cout << "chi2:\t"  << CHI2 << endl;
    const int ndof = hqcdp.NumberOfDegreesOfFreedom();
    cout << "Number of degrees of freedom (Ndof):\t" << ndof << endl;
    cout << "chi2/Ndof:\t" << CHI2/ndof << endl;

    return 0;
}