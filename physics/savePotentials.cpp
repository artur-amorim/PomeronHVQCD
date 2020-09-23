#include <iostream>
#include <exception>
#include "HolographicVQCD.h"
#include "schrodinger/schrodinger.h"
#include "methods/optimization/NelderMead.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    double sc, ksc, wsc, W0, w0, kU1, wU1;
    double VgIR, WIR, kIR, wIR, W1, k1, w1;
    double a1, a2, xf, tau0, Za, ca;
    string tag = "";
    if (argc < 22)
    {
        sc = 2.8328; ksc = 3.16697; wsc = 1.69263; W0 = 2.42885; w0 = 1.06381; kU1 = 1.32535;
        wU1 = -0.289843; VgIR = 1.80415; WIR = 1.13448; kIR = 1.76473; wIR = 2.92914; W1 = 0.234188; k1 = -0.307571;
        w1 = 3.03579; a1 = 0.541312; a2 = 1.54131; xf = 2.0/3; tau0 = 0.92323; Za = -0.0376777; ca = 23.3073;
    }
    else
    {
        sc = stod(argv[1]); ksc = stod(argv[2]); wsc = stod(argv[3]); W0 = stod(argv[4]); w0 = stod(argv[5]);
        kU1 = stod(argv[6]); wU1 = stod(argv[7]); VgIR = stod(argv[8]); WIR = stod(argv[9]); kIR = stod(argv[10]); 
        wIR = stod(argv[11]); W1 = stod(argv[12]); k1 = stod(argv[13]); w1 = stod(argv[14]); 
        a1 = stod(argv[15]); a2 = stod(argv[16]); xf = stod(argv[17]); tau0 = stod(argv[18]); Za = stod(argv[19]);
        ca = stod(argv[20]); tag = argv[21];
    }
    
    cout << "Computing potentials with parameter values" << endl;
    cout << "sc: " << sc << " ksc: " << ksc << " wsc: " << wsc << " W0: " << W0 << " w0: " << w0 << " kU1: " << kU1;
    cout << " wU1: " << wU1 << " VgIR: " << VgIR << " WIR: " << WIR << " kIR: " << kIR << " wIR: " << wIR << " W1: " << W1;
    cout << " k1: " << k1 << " w1: " << w1 <<  " a1: " << a1 << " a2: " << a2 << " xf: " << xf << " tau0: " << tau0;
    cout << " Za: " << Za << " ca: " << ca << endl;

    HVQCD hvqcd(sc, ksc, wsc, W0, w0, kU1, wU1, VgIR, WIR, kIR, wIR, W1, k1, w1, a1, a2, xf, tau0, Za, ca);
    hvqcd.solve(-10, 10);

    // Save background fields
    hvqcd.saveBackgroundFields("BackgroundFields_" + tag + ".txt");
    // Save background potentials profile
    hvqcd.savePotentials("BackgroundPotentials_" + tag + ".txt");
    // Save the schrodinger potentials
    saveSchrodingerPotentials(hvqcd, "SchrodingerPotentials_" + tag + ".txt");


    return 0;
}