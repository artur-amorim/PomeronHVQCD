#include <iostream>
#include <exception>
#include "HolographicVQCD.h"
#include "schrodinger/chebspec.h"
#include "schrodinger/schrodinger.h"
#include "methods/optimization/NelderMead.hpp"

using namespace std;

double J(const vector<double> X)
{
    double sc = 2.50485, ksc, wsc, W0, w0, kU1, wU1;
    double VgIR = 3.47852, WIR, kIR, wIR, W1, k1, w1;
    double a1, a2, xf = 2.0/3, tau0, Za = 133, ca = 0.26;
    ksc = X[0]; wsc = X[1]; W0 = X[2]; w0 = X[3]; kU1 = X[4];
    wU1 = X[5]; WIR = X[6]; kIR = X[7]; wIR = X[8]; W1 = X[9];
    k1 = X[10]; w1 = X[11]; a2 = X[12]; tau0 = X[13];
    a1 = a2*(4+W0*xf)/16.0;

    cout << "sc: " << sc << " ksc: " << ksc << " wsc: " << wsc << " W0: " << W0 << " w0: " << w0 << " kU1: " << kU1;
    cout << " wU1: " << wU1 << " VgIR: " << VgIR << " WIR: " << WIR << " kIR: " << kIR << " wIR: " << wIR << " W1: " << W1;
    cout << " k1: " << k1 << " w1: " << w1 << " a2: " << a2 << " tau0: " << tau0 << endl;
    
    HVQCD hvqcd(sc, ksc, wsc, W0, w0, kU1, wU1, VgIR, WIR, kIR, wIR, W1, k1, w1, a1, a2, xf, tau0, Za, ca);

    // Solve the background
    try
    {
        hvqcd.solve(-80, 20);
    }
    catch(...)
    {
        cout << "Bad background:" << endl;
        cout << "erms: " << 1e99 << endl;
        return 1e99;
    }

    // Compute the potentials
    vector<double> VVM, VAVM;
    try
    {
        VVM = computeVectorMesonPotential(hvqcd);
        VAVM = computeAxialVectorMesonNonSingletPotential(hvqcd, VVM);
    }
    catch(...)
    {
        cout << "Potentials not computed, check these values" << endl;
        cout << "erms: " << 1e99 << endl;
        return 1e99;
    }

    // Compute the masses
    vector<double> VMMasses, AVMMasses;
    try
    {
        vector<double> us = hvqcd.u();
        VMMasses = computeMasses(us, VVM, 6, "cheb");
        AVMMasses = computeMasses(us, VAVM, 5, "cheb");
    }
    catch(...)
    {
        cout << "Unable to compute spectrum." << endl;
        cout << "erms: " << 1e99 << endl;
        return 1e99;
    }
    double erms = 0;
    // In case we have negative m2 the mass container will be empty. In that case we fill it with zero entries
    // Contributions from the Meson ratios
    // Contributions from the ratios m_{\rho_n}/m_{\rho}
    if( VMMasses.size() == 0) VMMasses = vector<double>(6,0);
    for(int i = 0; i < Rrho_rho.size(); i++) erms += fabs((VMMasses[i+1]/VMMasses[0]-Rrho_rho[i])/Rrho_rho[i]);
    // Contributions from the ratios m_{\a1_n}/m_{\rho}
    if( AVMMasses.size() == 0) AVMMasses = vector<double>(5,0);
    for(int i = 0; i < Ra1_rho.size(); i++) erms += fabs((AVMMasses[i]/VMMasses[0]-Ra1_rho[i])/Ra1_rho[i]);
    // Singlet vector and axial vector meson sector
    for(int i = 0; i < Romega_rho.size(); i++) erms += fabs((VMMasses[i]/VMMasses[0]-Romega_rho[i])/Romega_rho[i]);
    //if( SingletAVMMasses.size() == 0) SingletAVMMasses = vector<double>(2,0);
    //for(int i = 0; i < Rf1_rho.size(); i++) erms += fabs((SingletAVMMasses[i]/VMMasses[0]-Rf1_rho[i])/Rf1_rho[i]);
    int nRatios = 13;
   
    erms = erms/nRatios;

    if (std::isnan(erms))
    {
        erms = 1e99;
        cout << "None erms" << endl;
        cout << "erms: " << erms << endl;
        return erms;
    }
    // Now we impose the constraint (12-x W0) kIR/VgIR/6>1
    double constr = (12-xf*W0)*kIR/(VgIR*6);
    erms += 0.1 * exp(- ( constr - 1) ); // The more the constraint is satisfied the smaller the penalty
    // We also impose the tachyo mass squared to be larger than 3.5
    double tmass2 = hvqcd.TachyonMassSquareIR();
    erms += 0.1 * exp( - (tmass2 - 3.5)) ;
    double mq = hvqcd.QuarkMass();
    cout << "(12 - xf) W0 kIR / (6 VgIR) = " << constr << " TachyonMassSquaredIR = " << tmass2 << endl;
    cout << "mq: " << mq << " erms: " << erms << endl;
    return erms;
}

int main(int argc, char ** argv)
{
    double sc = 2.50485, ksc, wsc, W0, w0, kU1, wU1;
    double VgIR = 3.47852, WIR, kIR, wIR, W1, k1, w1;
    double a2, xf = 2.0/3, tau0, Za = 133, ca = 0.26;
    if (argc < 15)
    {
        ksc = 3.0; wsc = 1.56; W0 = 2.5; w0 = 1.26; kU1 = 11./9; wU1 = 0.0;
        WIR = 0.9; kIR = 1.8; wIR = 5.0; W1 = 0.0; k1 = -0.23;
        w1 = 0.0; a2 = 1.0; tau0 = 1.;
    }
    else
    {
        ksc = stod(argv[1]); wsc = stod(argv[2]); W0 = stod(argv[3]); w0 = stod(argv[4]);
        kU1 = stod(argv[5]);
        wU1 = stod(argv[6]); WIR = stod(argv[7]); kIR = stod(argv[8]); wIR = stod(argv[9]);
        W1 = stod(argv[10]); k1 = stod(argv[11]); w1 = stod(argv[12]); a2 = stod(argv[13]);
        tau0 = stod(argv[14]);
    }
    
    cout << "Starting fit with values" << endl;
    cout << "sc: " << sc << " ksc: " << ksc << " wsc: " << wsc << " W0: " << W0 << " w0: " << w0 << " kU1: " << kU1;
    cout << " wU1: " << wU1 << " VgIR: " << VgIR << " WIR: " << WIR << " kIR: " << kIR << " wIR: " << wIR << " W1: " << W1;
    cout << " k1: " << k1 << " w1: " << w1 << " a2: " << a2 << " tau0: " << tau0 << endl;

    // Fit the model to the spectrum
    chebSetN(800);

    vector<double> x_guess = {ksc, wsc, W0, w0, kU1, wU1, WIR, kIR, wIR, W1, k1, w1, a2, tau0};

    vector<double> deltas = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    
    vector<double> xop = optimFunction(x_guess, J, deltas);
    
    // Show the optimum values found
    ksc = xop[0]; wsc = xop[1]; W0 = xop[2]; w0 = xop[3]; kU1 = xop[4];
    wU1 = xop[5]; WIR = xop[6]; kIR = xop[7]; wIR = xop[8];
    W1 = xop[9]; k1 = xop[10]; w1 = xop[11]; a2 = xop[12];
    tau0 = xop[13];
    double a1 = a2*(4+W0*xf)/16.0;

    HVQCD hvqcd(sc, ksc, wsc, W0, w0, kU1, wU1, VgIR, WIR, kIR, wIR, W1, k1, w1, a1, a2, xf, tau0, Za, ca);
    hvqcd.solve(-80,20);

    // Computing the mass ratios
    computeHVQCDRatios(hvqcd);

    double chi2 = J(xop);
    cout << "Best Chi2 found for ";
    cout << "sc: " << sc << " ksc: " << ksc << " wsc: " << wsc << " W0: " << W0 << " w0: " << w0 << " kU1: " << kU1;
    cout << " wU1: " << wU1 << " VgIR: " << VgIR << " WIR: " << WIR << " kIR: " << kIR << " wIR: " << wIR << " W1: " << W1;
    cout << " k1: " << k1 << " w1: " << w1 << " a1: " << a1 << " a2: " << a2 << " tau0: " << tau0 << endl;
    cout << "chi2: " << chi2 << endl;

    return 0;
}