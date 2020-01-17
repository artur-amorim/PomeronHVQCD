#include <iostream>
#include "HolographicVQCD.h"
#include "schrodinger/schrodinger.h"

using namespace std;

int main(int argc, char ** argv)
{
    double sc, ksc, wsc, W0, w0, kU1, wU1;
    double VgIR, WIR, kIR, wIR, W1, k1, w1;
    double xf, tau0;
    double za, ca;
    if (argc < 19)
    {
        cout << "Using default values" << endl;
        sc = 3; ksc = 3; wsc = 1.56; W0 = 2.5; w0 = 1.26; kU1 = 11.0/9; wU1 = 0.0;
        VgIR = 2.05; WIR = 0.9; kIR = 1.8; wIR = 5.0; W1 = 0.0; k1 = -0.23;
        w1 = 0.0; xf = 1.0; tau0 = 1.;
        za = 133; ca = 0.26;
    }
    else
    {
        sc = stod(argv[1]); ksc = stod(argv[2]); wsc = stod(argv[3]); W0 = stod(argv[4]); w0 = stod(argv[5]);
        kU1 = stod(argv[6]); wU1 = stod(argv[7]); VgIR = stod(argv[8]); WIR = stod(argv[9]); kIR = stod(argv[10]);
        wIR = stod(argv[11]); W1 = stod(argv[12]); k1 = stod(argv[13]); w1 = stod(argv[14]);
        xf = stod(argv[15]); tau0 = stod(argv[16]);
        za = stod(argv[17]); ca = stod(argv[18]);
    }
    // Create HVQCD object with default values
    HVQCD hvqcd(sc, ksc, wsc, W0, w0, kU1, wU1, VgIR, WIR, kIR, wIR, W1, k1, w1, xf, tau0, za, ca, false, true, true, true, true);
    // Solve the background
    hvqcd.solve();
    hvqcd.printQuarkMass();
    // Print the ratio values
    chebSetN(800);
    hvqcd.computeSpectrum();
    hvqcd.showRatioValues();

    double tmass2 = hvqcd.TachyonMassSquareIR();

    cout << "Tau Mass Squered IR: " << tmass2 << endl;

    return 0;
}