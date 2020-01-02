#include <exception>
#include "HolographicVQCD.h"
#include "schrodinger/schrodinger.h"

using namespace std;

int main(int argc, char ** argv)
{
    double sc, ksc, wsc, W0, w0, wU1;
    double VgIR, WIR, kIR, wIR, W1, k1, w1;
    double xf, tau0;
    double Za, ca;
    bool add_quark, add_tensor_glueball, add_scalars, add_singlet_vector, add_singlet_axial;
    if (argc < 22)
    {
        sc = 3; ksc = 3; wsc = 1.56; W0 = 2.5; w0 = 1.26; wU1 = 0.0;
        VgIR = 2.05; WIR = 0.9; kIR = 1.8; wIR = 5.0; W1 = 0.0; k1 = -0.23;
        w1 = 0.0; xf = 1.0; tau0 = 1.; 
        Za = 133; ca = 0.26;
        add_quark = false; add_tensor_glueball = true;
        add_scalars = true; add_singlet_vector = false; add_singlet_axial = false;
    }
    else
    {
        sc = stod(argv[1]); ksc = stod(argv[2]); wsc = stod(argv[3]); W0 = stod(argv[4]); w0 = stod(argv[5]);
        wU1 = stod(argv[6]); VgIR = stod(argv[7]); WIR = stod(argv[8]); kIR = stod(argv[9]); wIR = stod(argv[10]);
        W1 = stod(argv[11]); k1 = stod(argv[12]); w1 = stod(argv[13]); xf = 2.0/3; tau0 = stod(argv[14]);
        Za = stod(argv[15]); ca = stod(argv[16]);
        int i1 = stoi(argv[17]), i2 = stoi(argv[18]), i3 = stoi(argv[19]), i4 = stoi(argv[20]), i5 = stoi(argv[21]);
        if (i1 == 0 || i1 == 1) add_quark = bool(i1);
        else throw domain_error("add_quark is neither 0 or 1.");
        if (i2 == 0 || i2 == 1) add_tensor_glueball = bool(i2);
        else throw domain_error("add_tensor_glueball is neither 0 or 1.");
        if (i3 == 0 || i3 == 1) add_scalars = bool(i3);
        else throw domain_error("add_scalars is neither 0 or 1.");
        if (i4 == 0 || i4 == 1) add_singlet_vector = bool(i4);
        else throw domain_error("add_singlet_vector is neither 0 or 1.");
        if (i5 == 0 || i5 == 1) add_singlet_axial = bool(i5);
        else throw domain_error("add_singlet_axial is neither 0 or 1.");
    }

    cout << "Computing mass ratios for" << endl; 
    cout << "sc: " << sc << " ksc: " << ksc << " wsc: " << wsc << " W0: " << W0 << " w0: " << w0;
    cout << " wU1: " << wU1 << " VgIR: " << VgIR << " WIR: " << WIR << " kIR: " << kIR << " wIR: " << wIR << " W1: " << W1;
    cout << " k1: " << k1 << " w1: " << w1 << " tau0: " << tau0 << " Za: " << Za << " ca: " << ca << endl;
    // Create HVQCD object
    HVQCD hvqcd(sc, ksc, wsc, W0, w0, 11./9, wU1, VgIR, WIR, kIR,
                wIR, W1, k1, w1, xf, tau0, Za, ca, add_quark, add_tensor_glueball, add_scalars, add_singlet_vector, add_singlet_axial);
    // Fit the model to the spectrum
    chebSetN(800);
    hvqcd.showRatioValues();

    // Save background field values
    hvqcd.saveBackgroundFields("plots/HVQCD/BackgroundFields.txt");
    // Now let's save the values of w, k, Vg, W and so on
    hvqcd.savePotentials("plots/HVQCD/Potentials.txt");

    return 0;
}