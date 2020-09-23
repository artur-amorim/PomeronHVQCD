#include <iostream>
#include <vector>
#include <functional>
#include "HolographicVQCD.h"
#include "F2.h"
#include "FL.h"
#include "GluonKernel.h"
#include "MesonKernel.h"
#include "HQCDP.h"
#include "schrodinger/schrodinger.h"
#include "schrodinger/chebspec.h"
#include "methods/optimization/NelderMead.hpp"


using namespace std;

int main(int argc, char ** argv)
{
    double invls, ag, bg, cg, dg, eg, fg;
    double am, bm, cm, dm, em, fm;
    double g1g, g2g, g1m;
    if (argc < 17)
    {
        invls = -0.00329973; ag = -4.42448; bg = -11.9668; cg = 2.00837; dg = 2.75436; eg = -0.229314; fg = 4.79135;
        am = 4.5345; bm = -4.46109; cm = 10.2159; dm = 0.444645; em = 6.65785; fm = 1.16794;
        g1g = 0; g2g = 0; g1m = 0;  
    }
    else
    {
        invls = stod(argv[1]); ag = stod(argv[2]); bg = stod(argv[3]); cg = stod(argv[4]); dg = stod(argv[5]);
        eg = stod(argv[6]); fg = stod(argv[7]);
        am = stod(argv[8]); bm = stod(argv[9]); cm = stod(argv[10]); dm = stod(argv[11]);
        em = stod(argv[12]); fm = stod(argv[13]);
        g1g = stod(argv[14]); g2g = stod(argv[15]); g1m = stod(argv[16]);
    }

    cout << "Starting the search of gs with" << endl;
    cout << "invls: " << invls << " ag: " << ag << " bg: " << bg << " cg: " << cg << " dg: " << dg;
    cout << " eg: " << eg <<  " fg: " << fg << endl;
    cout << "am: " << am << " bm: " << bm << " cm: " << cm << " dm: " << dm;
    cout << " em: " << em <<  " fm: " << fm << endl;
    cout << "g1g: " << g1g << " g2g: " << g2g << " g1m: " << g1m << endl;

    double mq = hvqcd().QuarkMass();
    F2 f2("expdata/DIS/F2_data.txt");
    FL fl("expdata/DIS/FL_data.txt");

    // Setup Gluon Kernel and GNs vector
    GluonKernel gluon(2, {invls, ag, bg, cg, dg, eg, fg});
    MesonKernel meson(1, {invls, am, bm, cm, dm, em, fm});
    vector<double> GNs = {g1g, g2g, g1m};

    // Setup HQCDP object
    HQCDP hqcdp;
    hqcdp.addProcessObservable(f2);
    hqcdp.addProcessObservable(fl);
    hqcdp.addKernel(gluon);
    hqcdp.addKernel(meson);
    hqcdp.setGNs(GNs);
    chebSetN(400);
    hqcdp.computeSpectrum();

    // Function that we want to minimize
    function<double(vector<double>)> func = [&hqcdp] (const vector<double> &pars)
    {
        hqcdp.setGNs(pars);
        cout << "g1g : " << pars[0] << " g2g: " << pars[1] << " g1m: " << pars[2] << endl;
        double chi2 = hqcdp.chi2();
        cout << "chi2: " << chi2 << endl;
        return chi2;
    };

    GNs = optimFunction(GNs, func, 10, 1e-12);

    cout << "Best gns found are:" << endl;
    cout << "g1g:\t" << GNs[0] << endl;
    cout << "g2g:\t" << GNs[1] << endl;
    cout << "g1m:\t" << GNs[2] << endl;

    return 0;
}