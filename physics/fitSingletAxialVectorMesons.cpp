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
        HVQCD hvqcd(1.76662, 2.65341, 0.0258923, 2.74747, 1.39882, 11./9, 6.25665, 0.954371, 2.73599, 0.294752, 12.6054,
                1.2582, 5.07657, 0.75851, 2./3, 3.03924, x[0], x[1], add_quark, add_tensor_glueball, add_scalars, add_singlet_vector, add_singlet_axial);
        double j = hvqcd.J();
        return j;
    }
};

int main(int argc, char ** argv)
{
    double sc = 1.76662, ksc = 2.65341, wsc = 0.0258923, W0 = 2.74747, w0 = 1.39882, wU1 = 6.25665;
    double VgIR = 0.954371, WIR = 2.73599, kIR = 0.294752, wIR = 12.6054, W1 = 1.2582, k1 = 5.07657, w1 = 0.75851;
    double kU1 = 11./9, xf = 2.0/3, tau0 = 3.03924, Za = 3.07123e-07, ca = 8.67959e-06;
    double deltaZa = 1e-07, deltaca = 1e-06;
    bool add_quark = false, add_tensor_glueball = true, add_scalars = false;
    bool add_singlet_vector = true, add_singlet_axial = true;
    if (argc < 5)
    {
        cout << "Starting fit with default values:" << endl;
        cout << "Za: " << Za << " ca: " << ca << endl;
        cout << "deltaZa: " << deltaZa << " deltaca: " << deltaca << endl;
    }
    else
    {
        Za = stod(argv[1]); ca = stod(argv[2]);
        deltaZa = stod(argv[3]); ca = stod(argv[4]);
    }

    // Fit the model to the spectrum
    chebSetN(800);

    vector<double> x_guess = {Za, ca};
    vector<double> deltas = {deltaZa, deltaca};

    optmFunc func(add_quark, add_tensor_glueball, add_scalars, add_singlet_vector, add_singlet_axial);
    
    Fit bestFit = optimFunction(x_guess, func, deltas);
    vector<double> xop = bestFit.bestFitParsUncert;

    // Show the optimum values found
    Za = xop[0]; ca = xop[1];

    HVQCD hvqcd(sc, ksc, wsc, W0, w0, kU1, wU1, VgIR, WIR, kIR, wIR, W1, k1, w1, xf, tau0, Za, ca,
            add_quark, add_tensor_glueball, add_scalars, add_singlet_vector, add_singlet_axial);

    double chi2 = hvqcd.J();
    cout << "Best Chi2 found for ";
    cout << "sc: " << sc << " ksc: " << ksc << " wsc: " << wsc << " W0: " << W0 << " w0: " << w0;
    cout << " wU1: " << wU1 << " VgIR: " << VgIR << " WIR: " << WIR << " kIR: " << kIR << " wIR: " << wIR << " W1: " << W1;
    cout << " k1: " << k1 << " w1: " << w1 << " tau0: " << tau0 <<  " Za: " << Za << " ca: " << ca << endl;
    cout << "chi2: " << chi2 << endl;

    // Computing the mass ratios
    hvqcd.showRatioValues();

    return 0;
}