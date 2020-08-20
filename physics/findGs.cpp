#include <iostream>
#include <vector>
#include <functional>
#include "HolographicVQCD.h"
#include "F2.h"
#include "FL.h"
#include "GluonKernel.h"
#include "HQCDP.h"
#include "schrodinger/schrodinger.h"
#include "schrodinger/chebspec.h"
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
        g1 = 0; g2 = 0; g3 = 0; g4 = 0;  
    }
    else
    {
        invls = stod(argv[1]); a = stod(argv[2]); b = stod(argv[3]); c = stod(argv[4]); d= stod(argv[5]);
        e = stod(argv[6]); f = stod(argv[7]);
        g1 = stod(argv[8]); g2 = stod(argv[9]); g3 = stod(argv[10]); g4 = stod(argv[11]);
    }

    cout << "Starting the search of gs with" << endl;
    cout << "invls: " << invls << " a: " << a << " b: " << b << " c: " << c << " d: " << d;
    cout << " e: " << e <<  " f: " << f << endl;
    cout << " g1: " << g1 << " g2: " << g2 << " g3: " << g3 << " g4: " << g4 << endl;

    double mq = hvqcd().QuarkMass();
    F2 f2("expdata/F2_data.txt");
    FL fl("expdata/FL_data.txt");

    // Setup Gluon Kernel and GNs vector
    vector<double> gluon_pars = {invls, a, b, c, d, e, f};
    GluonKernel gluon(4, gluon_pars);
    vector<double> GNs = {g1, g2, g3, g4};

    // Setup HQCDP object
    HQCDP hqcdp;
    hqcdp.addProcessObservable(f2);
    hqcdp.addProcessObservable(fl);
    hqcdp.addKernel(gluon);
    hqcdp.setGNs(GNs);
    chebSetN(400);
    hqcdp.computeSpectrum();

    // Function that we want to minimize
    function<double(vector<double>)> func = [&hqcdp] (const vector<double> &pars)
    {
        hqcdp.setGNs(pars);
        cout << "g1 : " << pars[0] << " g2: " << pars[1] << " g3: " << pars[2] << " g4: " << pars[3] << endl;
        double chi2 = hqcdp.chi2();
        cout << "chi2: " << chi2 << endl;
        return chi2;
    };

    GNs = optimFunction(GNs, func, 10, 1e-3);

    cout << "Best gns found are:" << endl;
    cout << "g0:\t" << GNs[0] << endl;
    cout << "g1:\t" << GNs[1] << endl;
    cout << "g2:\t" << GNs[2] << endl;
    cout << "g3:\t" << GNs[3] << endl;

    return 0;
}