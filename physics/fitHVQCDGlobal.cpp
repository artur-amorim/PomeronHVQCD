#include <iostream>
#include <exception>
#include "HolographicVQCD.h"
#include "schrodinger/schrodinger.h"
#include "methods/optimization/NelderMead.hpp"

using namespace std;

double J(const vector<double> X)
{
    double sc, ksc, wsc, W0, w0, kU1, wU1;
    double VgIR, WIR, kIR, wIR, W1, k1, w1;
    double xf = 2.0/3, tau0, Za, ca;
    sc = X[0]; ksc = X[1]; wsc = X[2]; W0 = X[3]; w0 = X[4];
    kU1 = X[5];
    wU1 = X[6]; VgIR = X[7]; WIR = X[8]; kIR = X[9]; wIR = X[10];
    W1 = X[11]; k1 = X[12]; w1 = X[13]; tau0 = X[14];
    Za = X[15]; ca = X[16];

    cout << "sc: " << sc << " ksc: " << ksc << " wsc: " << wsc << " W0: " << W0 << " w0: " << w0 << " kU1: " << kU1;
    cout << " wU1: " << wU1 << " VgIR: " << VgIR << " WIR: " << WIR << " kIR: " << kIR << " wIR: " << wIR << " W1: " << W1;
    cout << " k1: " << k1 << " w1: " << w1 << " tau0: " << tau0 <<  " Za: " << Za << " ca: " << ca << endl;
    
    HVQCD hvqcd(sc, ksc, wsc, W0, w0, kU1, wU1, VgIR, WIR, kIR, wIR, W1, k1, w1, xf, tau0, Za, ca);

    // Solve the background
    try
    {
        hvqcd.solve();
    }
    catch(...)
    {
        cout << "Bad background:" << endl;
        cout << "erms: " << 1e99 << endl;
        return 1e99;
    }

    // Compute the potentials
    vector<double> V2G, VVM, VAVM, VPSM, VSM, VSingletAVM;
    try
    {
        V2G = computeV2GPotential(hvqcd);
        VVM = computeVectorMesonPotential(hvqcd);
        VAVM = computeAxialVectorMesonNonSingletPotential(hvqcd, VVM);
        VPSM = computePseudoScalarMesonPotential(hvqcd);
        VSM = computeScalarMesonPotential(hvqcd);
        VSingletAVM = computeAxialVectorMesonSingletPotential(hvqcd, VAVM);
    }
    catch(...)
    {
        cout << "Potentials not computed, check these values" << endl;
        cout << "erms: " << 1e99 << endl;
        return 1e99;
    }

    // Compute the masses
    vector<double> TGMasses, VMMasses, AVMMasses, PSMMasses, SMMasses, SingletAVMMasses;
    try
    {
        vector<double> zs = hvqcd.z(), us = hvqcd.u();
        TGMasses = computeMasses(zs, V2G, 1, "cheb");
        VMMasses = computeMasses(us, VVM, 6, "cheb");
        AVMMasses = computeMasses(us, VAVM, 5, "cheb");
        PSMMasses = computeMasses(us, VPSM, 5, "cheb");
        SMMasses = computeMasses(us,VSM, 2, "cheb");
        SingletAVMMasses = computeMasses(us, VSingletAVM, 2, "cheb");
        if(SingletAVMMasses.size() == 0) throw(runtime_error("Only negative values for Singlet AVM masses"));
    }
    catch(...)
    {
        cout << "Unable to compute spectrum. Check these values" << endl;
        cout << "erms: " << 1e99 << endl;
        return 1e99;
    }
    double erms = 0;
    // Chi2 contribution due to the tensor glueball ratio
    erms += pow((TGMasses[0]/VMMasses[0] - RTG_rho[0])/RTG_rho[0], 2);
    // Contributions from the Meson ratios
    // Contributions from the ratios m_{\rho_n}/m_{\rho}
    for(int i = 0; i < Rrho_rho.size(); i++) erms += pow((VMMasses[i+1]/VMMasses[0]-Rrho_rho[i])/Rrho_rho[i],2);
    // Contributions from the ratios m_{\a1_n}/m_{\rho}
    for(int i = 0; i < Ra1_rho.size(); i++) erms += pow((AVMMasses[i]/VMMasses[0]-Ra1_rho[i])/Ra1_rho[i],2);
    // Contributions from the ratios m_{\pi_n} / m_{\rho}
    for(int i = 0; i < Rpi_rho.size(); i++) erms += pow((PSMMasses[i]/VMMasses[0]-Rpi_rho[i])/Rpi_rho[i],2);
    // Contributions from the ratios m_{a0_n} / m_{\rho}
    for(int i = 0; i < Ra0_rho.size(); i++) erms += pow((SMMasses[i]/VMMasses[0]- Ra0_rho[i])/Ra0_rho[i], 2);
    // Singlet vector and axial vector meson sector
    for(int i = 0; i < Romega_rho.size(); i++) erms += pow((VMMasses[i]/VMMasses[0]-Romega_rho[i])/Romega_rho[i],2);
    for(int i = 0; i < Rf1_rho.size(); i++) erms += pow((SingletAVMMasses[i]/VMMasses[0]-Rf1_rho[i])/Rf1_rho[i],2);
    int nRatios = 23;
   
    erms = sqrt(erms)/nRatios;

    if (std::isnan(erms))
    {
        erms = 1e99;
        cout << "None erms" << endl;
        cout << "erms: " << erms << endl;
        return erms;
    }
    // Now we impose the constraint (12-x W0) kIR/VgIR/6>1
    double constr = (12-xf*W0)*kIR/(VgIR*6);
    erms += exp(- ( constr - 1) ); // The more the constraint is satisfied the smaller the penalty
    // We also impose the tachyo mass squared to be larger than 3.5
    double tmass2 = hvqcd.TachyonMassSquareIR();
    erms += exp( - (tmass2 - 3.5)) ;
    double mq = hvqcd.QuarkMass();
    cout << "(12 - xf) W0 kIR / (6 VgIR) = " << constr << " TachyonMassSquaredIR = " << tmass2 << endl;
    cout << "mq: " << mq << " erms: " << erms << endl;
    return erms;
}

int main(int argc, char ** argv)
{
    double sc, ksc, wsc, W0, w0, kU1, wU1;
    double VgIR, WIR, kIR, wIR, W1, k1, w1;
    double xf = 2.0/3, tau0, Za, ca;
    if (argc < 18)
    {
        sc = 3.0; ksc = 3.0; wsc = 1.56; W0 = 2.5; w0 = 1.26; kU1 = 11./9; wU1 = 0.0;
        VgIR = 2.05; WIR = 0.9; kIR = 1.8; wIR = 5.0; W1 = 0.0; k1 = -0.23;
        w1 = 0.0; tau0 = 1.; Za = 20; ca = 0.26;
    }
    else
    {
        sc = stod(argv[1]); ksc = stod(argv[2]); wsc = stod(argv[3]); W0 = stod(argv[4]); w0 = stod(argv[5]);
        kU1 = stod(argv[6]);
        wU1 = stod(argv[7]); VgIR = stod(argv[8]); WIR = stod(argv[9]); kIR = stod(argv[10]); wIR = stod(argv[11]);
        W1 = stod(argv[12]); k1 = stod(argv[13]); w1 = stod(argv[14]); tau0 = stod(argv[15]);
        Za = stod(argv[16]); ca = stod(argv[17]);
    }
    
    cout << "Starting fit with values" << endl;
    cout << "sc: " << sc << " ksc: " << ksc << " wsc: " << wsc << " W0: " << W0 << " w0: " << w0 << " kU1: " << kU1;
    cout << " wU1: " << wU1 << " VgIR: " << VgIR << " WIR: " << WIR << " kIR: " << kIR << " wIR: " << wIR << " W1: " << W1;
    cout << " k1: " << k1 << " w1: " << w1 << " tau0: " << tau0 <<  " Za: " << Za << " ca: " << ca << endl;

    // Fit the model to the spectrum
    chebSetN(800);

    vector<double> x_guess = {sc, ksc, wsc, W0, w0, kU1, wU1, VgIR, WIR, kIR, wIR, W1, k1, w1, tau0, Za, ca};

    vector<double> deltas = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 10.0, 0.01};
    
    vector<double> xop = optimFunction(x_guess, J, deltas);
    
    // Show the optimum values found
    sc = xop[0]; ksc = xop[1]; wsc = xop[2]; W0 = xop[3]; w0 = xop[4]; kU1 = xop[5];
    wU1 = xop[6]; VgIR = xop[7]; WIR = xop[8]; kIR = xop[9]; wIR = xop[10];
    W1 = xop[11]; k1 = xop[12]; w1 = xop[13]; xf = 2./3; tau0 = xop[14];
    Za = xop[15]; ca = xop[16];

    HVQCD hvqcd(sc, ksc, wsc, W0, w0, kU1, wU1, VgIR, WIR, kIR, wIR, W1, k1, w1, xf, tau0, Za, ca);
    hvqcd.solve();

    // Computing the mass ratios
    computeHVQCDRatios(hvqcd);

    double chi2 = J(xop);
    cout << "Best Chi2 found for ";
    cout << "sc: " << sc << " ksc: " << ksc << " wsc: " << wsc << " W0: " << W0 << " w0: " << w0 << " kU1: " << kU1;
    cout << " wU1: " << wU1 << " VgIR: " << VgIR << " WIR: " << WIR << " kIR: " << kIR << " wIR: " << wIR << " W1: " << W1;
    cout << " k1: " << k1 << " w1: " << w1 << " tau0: " << tau0 <<  " Za: " << Za << " ca: " << ca << endl;
    cout << "chi2: " << chi2 << endl;

    return 0;
}