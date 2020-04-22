#include <iostream>
#include <vector>
#include "HolographicVQCD.h"
#include "F2.h"
#include "FL.h"
#include "GluonKernel.h"
#include "HQCDP.h"
#include "schrodinger/schrodinger.h"


using namespace std;

int main(int argc, char ** argv)
{
    double invls, a, b, c, d, e, f;
    double g1, g2, g3, g4;
    if (argc < 12)
    {
        invls = 1.0/0.153; a = 0; b = 0; c = 0; d = 0; e = 0; f = 0;
        g1 = 0; g2 = 0; g3 = 0; g4 = 0;  
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
    cout << "mq: " << mq << endl;
    F2 f2(false, "expdata/F2_data.txt");
    FL fl(false, "expdata/FL_data.txt");

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
    
    // Compute the spectrum to check Reggeon properties
    chebSetN(1000);

    vector<double> guess = {invls, a, b, c, d, e, f, g1, g2, g3, g4};
    hqcdp.fit(guess, 0.5);

    return 0;
}