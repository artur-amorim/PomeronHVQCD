#include <iostream>
#include "HolographicVQCD.h"
#include "schrodinger/schrodinger.h"

using namespace std;

int main(int argc, char ** argv)
{
    double sc, ksc, wsc, W0, w0, kU1, wU1;
    double VgIR, WIR, kIR, wIR, W1, k1, w1;
    double a1, a2;
    double xf, tau0;
    double za, ca;
    if (argc < 19)
    {
        cout << "Using default values" << endl;
        sc = 3; ksc = 3; wsc = 1.56; W0 = 2.5; w0 = 1.26; kU1 = 11.0/9; wU1 = 0.0;
        VgIR = 2.05; WIR = 0.9; kIR = 1.8; wIR = 5.0; W1 = 0.0; k1 = -0.23;
        w1 = 0.0; a1 = 0.0; a2 = 1.0; xf = 1.0; tau0 = 1.; za = 133; ca = 0.26;
    }
    else
    {
        sc = stod(argv[1]); ksc = stod(argv[2]); wsc = stod(argv[3]); W0 = stod(argv[4]); w0 = stod(argv[5]);
        kU1 = stod(argv[6]); wU1 = stod(argv[7]); VgIR = stod(argv[8]); WIR = stod(argv[9]); kIR = stod(argv[10]);
        wIR = stod(argv[11]); W1 = stod(argv[12]); k1 = stod(argv[13]); w1 = stod(argv[14]); a1 = stod(argv[15]);
        a2 = stod(argv[16]); xf = stod(argv[17]); tau0 = stod(argv[18]); za = stod(argv[19]); ca = stod(argv[20]);
    }
    cout << "sc: " << sc << " ksc: " << ksc << " wsc: " << wsc << " W0: " << W0 << " w0: " << w0 << endl;
    cout << "kU1: " << kU1 << " wU1: " << wU1 << " VgIR: " << VgIR << " WIR: " << WIR << " kIR: " << kIR << endl;
    cout << "wIR: " << wIR << " W1: " << W1 << " k1: " << k1 << " w1: " << w1 << " a1: " << a1 << " a2: " << a2 << endl;
    cout << "x: " << xf << " tau0: " << tau0 << " za: " << za << " ca: " << ca << endl;
    // Create HVQCD object with default values
    HVQCD hvqcd(sc, ksc, wsc, W0, w0, kU1, wU1, VgIR, WIR, kIR, wIR, W1, k1, w1, a1, a2, xf, tau0, za, ca);
    // Solve the background
    hvqcd.solve();
    
    double mq = hvqcd.QuarkMass();
    double tmass2 = hvqcd.TachyonMassSquareIR();

    cout << "Quark Mass: " << mq << endl;
    cout << "Tau Mass Squered IR: " << tmass2 << endl;

    hvqcd.saveBackgroundFields("BackgroundFields.txt");
    saveSchrodingerPotentials(hvqcd);

    chebSetN(800);
    computeHVQCDSpectrum(hvqcd);
    computeHVQCDRatios(hvqcd);

    return 0;
}