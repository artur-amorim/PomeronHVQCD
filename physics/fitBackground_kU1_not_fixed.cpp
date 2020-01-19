#include <iostream>
#include <exception>
#include "HolographicVQCD.h"
#include "schrodinger/schrodinger.h"
#include "methods/optimization/NelderMead.hpp"

using namespace std;

struct optmFunc
{
    bool add_quark, add_tensor_glueball, add_scalars;
    bool add_singlet_vector, add_singlet_axial;
    optmFunc(): add_quark(false), add_tensor_glueball(true), add_scalars(false),
                add_singlet_vector(false), add_singlet_axial(false) {}
    optmFunc(const bool quark, const bool tg, const bool scal, const bool sing_vec, const bool sing_ax):
            add_quark(quark), add_tensor_glueball(tg), add_scalars(scal),
            add_singlet_vector(sing_vec), add_singlet_axial(sing_ax) {}
    double operator()(const std::vector<double> &x)
    {
        HVQCD hvqcd(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10],
                x[11], x[12], x[13], 2./3, x[14], 0, 0, add_quark, add_tensor_glueball, add_scalars, add_singlet_vector, add_singlet_axial);
        double j = hvqcd.J();
        return j;
    }
};

int main(int argc, char ** argv)
{
    double sc, ksc, wsc, W0, w0, kU1, wU1;
    double VgIR, WIR, kIR, wIR, W1, k1, w1;
    double xf = 2.0/3, tau0, Za, ca;
    bool add_quark = false, add_tensor_glueball = true, add_scalars = false;
    bool add_singlet_vector = true, add_singlet_axial = false;
    if (argc < 18)
    {
        cout << "Starting fit with default values" << endl;
        sc = 3; ksc = 3; wsc = 1.56; W0 = 2.5; w0 = 1.26; kU1 = 11./9; wU1 = 0.0;
        VgIR = 2.05; WIR = 0.9; kIR = 1.8; wIR = 5.0; W1 = 0.0; k1 = -0.23;
        w1 = 0.0; tau0 = 1.; Za = 133; ca = 0.26;
    }
    else
    {
        sc = stod(argv[1]); ksc = stod(argv[2]); wsc = stod(argv[3]); W0 = stod(argv[4]); w0 = stod(argv[5]);
        kU1 = stod(argv[6]);
        wU1 = stod(argv[7]); VgIR = stod(argv[8]); WIR = stod(argv[9]); kIR = stod(argv[10]); wIR = stod(argv[11]);
        W1 = stod(argv[12]); k1 = stod(argv[13]); w1 = stod(argv[14]); tau0 = stod(argv[15]);
        Za = stod(argv[16]); ca = stod(argv[17]);
    }

    // Fit the model to the spectrum
    chebSetN(800);

    vector<double> x_guess = {sc, ksc, wsc, W0, w0, kU1, wU1, VgIR, WIR, kIR, wIR, W1, k1, w1, tau0};

    optmFunc func(add_quark, add_tensor_glueball, add_scalars, add_singlet_vector, add_singlet_axial);

    vector<double> deltas = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    
    Fit bestFit = optimFunction(x_guess, func, deltas);
    vector<double> xop = bestFit.bestFitPars;
    
    // Show the optimum values found
    sc = xop[0]; ksc = xop[1]; wsc = xop[2]; W0 = xop[3]; w0 = xop[4]; kU1 = xop[5];
    wU1 = xop[6]; VgIR = xop[7]; WIR = xop[8]; kIR = xop[9]; wIR = xop[10];
    W1 = xop[11]; k1 = xop[12]; w1 = xop[13]; xf = 2./3; tau0 = xop[14];
    //Za = xop[15]; ca = xop[16];

    HVQCD hvqcd(sc, ksc, wsc, W0, w0, kU1, wU1, VgIR, WIR, kIR, wIR, W1, k1, w1, xf, tau0, Za, ca,
            add_quark, add_tensor_glueball, add_scalars, add_singlet_vector, add_singlet_axial);

    double chi2 = hvqcd.J();
    cout << "Best Chi2 found for ";
    cout << "sc: " << sc << " ksc: " << ksc << " wsc: " << wsc << " W0: " << W0 << " w0: " << w0 << " kU1: " << kU1;
    cout << " wU1: " << wU1 << " VgIR: " << VgIR << " WIR: " << WIR << " kIR: " << kIR << " wIR: " << wIR << " W1: " << W1;
    cout << " k1: " << k1 << " w1: " << w1 << " tau0: " << tau0 <<  " Za: " << Za << " ca: " << ca << endl;
    cout << "chi2: " << chi2 << endl;

    // Computing the mass ratios
    hvqcd.showRatioValues();

    return 0;
}