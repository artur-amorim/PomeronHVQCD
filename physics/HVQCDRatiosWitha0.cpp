#include <iostream>
#include <exception>
#include "HolographicVQCD.h"
#include "schrodinger/chebspec.h"
#include "schrodinger/schrodinger.h"

using namespace std;

int main(int argc, char ** argv)
{
    double sc, ksc, wsc, W0, w0, kU1, wU1;
    double VgIR, WIR, kIR, wIR, W1, k1, w1;
    double a1, a2, xf = 2.0/3, tau0, Za, ca;
    if (argc < 18)
    {
        sc = 3.0; ksc = 3.0; wsc = 1.56; W0 = 2.5; w0 = 1.26; kU1 = 11./9; wU1 = 0.0;
        VgIR = 2.05; WIR = 0.9; kIR = 1.8; wIR = 5.0; W1 = 0.0; k1 = -0.23;
        w1 = 0.0; a1 = 0; a2 = 1; tau0 = 1.; Za = 133; ca = 0.26;
    }
    else
    {
        sc = stod(argv[1]); ksc = stod(argv[2]); wsc = stod(argv[3]); W0 = stod(argv[4]); w0 = stod(argv[5]);
        kU1 = stod(argv[6]);
        wU1 = stod(argv[7]); VgIR = stod(argv[8]); WIR = stod(argv[9]); kIR = stod(argv[10]); wIR = stod(argv[11]);
        W1 = stod(argv[12]); k1 = stod(argv[13]); w1 = stod(argv[14]); a1 = stod(argv[15]); a2 = stod(argv[16]); tau0 = stod(argv[17]);
        Za = stod(argv[18]); ca = stod(argv[19]);
    }
    
    cout << "Computing ratios with parameter values" << endl;
    cout << "sc: " << sc << " ksc: " << ksc << " wsc: " << wsc << " W0: " << W0 << " w0: " << w0 << " kU1: " << kU1;
    cout << " wU1: " << wU1 << " VgIR: " << VgIR << " WIR: " << WIR << " kIR: " << kIR << " wIR: " << wIR << " W1: " << W1;
    cout << " k1: " << k1 << " w1: " << w1 << " a1: " << a1 << " a2: " << a2 << " tau0: " << tau0 << " Za: " << Za << " ca: " << ca << endl;

    // Fit the model to the spectrum
    chebSetN(800);

    HVQCD hvqcd(sc, ksc, wsc, W0, w0, kU1, wU1, VgIR, WIR, kIR, wIR, W1, k1, w1, a1, a2, xf, tau0, Za, ca);
    hvqcd.solve(-80, 20);

    const vector<double> zs = hvqcd.z(), us = hvqcd.u();

    // Compute the Schrodinger potential of the fluctuations
    std::vector<double> V2G, VVM, VAVM, VSM, VSingletAVM;
    V2G = computeV2GPotential(hvqcd);
    VVM = computeVectorMesonPotential(hvqcd);
    VAVM = computeAxialVectorMesonNonSingletPotential(hvqcd, VVM);
    VSM = computeScalarMesonPotential(hvqcd);
    VSingletAVM = computeAxialVectorMesonSingletPotential(hvqcd, VAVM);
    // Compute the masses
    std::vector<double> TGMasses, VMMasses, AVMMasses, PSMMasses, SMMasses, SingletAVMMasses;
    TGMasses = computeMasses(zs, V2G, 1);
    VMMasses = computeMasses(us, VVM, 7);
    AVMMasses = computeMasses(us, VAVM, 4);
    PSMMasses = computePseudoScalarMasses(hvqcd, 5);
    SMMasses = computeMasses(us,VSM, 3);
    SingletAVMMasses = computeMasses(us, VSingletAVM, 4);
    // Compute the predicted ratios
    std::vector<double> RTGPred, RrhoPred, Ra1Pred, RpiPred, Ra0Pred, RomegaPred, Rf1Pred;
    RTGPred = {TGMasses[0]/ VMMasses[0]};
    RrhoPred = {VMMasses[1]/VMMasses[0], VMMasses[2]/VMMasses[0], VMMasses[3]/VMMasses[0], VMMasses[4]/VMMasses[0]};
    Ra1Pred = {AVMMasses[0]/VMMasses[0], AVMMasses[1]/VMMasses[0], AVMMasses[2]/VMMasses[0], AVMMasses[3]/VMMasses[0]};
    RpiPred = {PSMMasses[0]/VMMasses[0], PSMMasses[1]/VMMasses[0], PSMMasses[2]/VMMasses[0], PSMMasses[3]/VMMasses[0], PSMMasses[4]/VMMasses[0]};
    Ra0Pred = {SMMasses[0] / VMMasses[0], SMMasses[1] / VMMasses[0], SMMasses[2] / VMMasses[0]};
    RomegaPred = {VMMasses[0]/VMMasses[0], VMMasses[1]/VMMasses[0], VMMasses[2]/VMMasses[0], VMMasses[3]/VMMasses[0], VMMasses[4]/VMMasses[0], VMMasses[5]/VMMasses[0], VMMasses[6]/VMMasses[0]};
    Rf1Pred = {SingletAVMMasses[0]/VMMasses[0], SingletAVMMasses[1]/VMMasses[0], SingletAVMMasses[2]/VMMasses[0], SingletAVMMasses[3]/VMMasses[0]};
    double sum_errors = 0;
    // Compare the predicted ratios with the known ones
    const std::vector<double> new_Ra0_rho = {1.2640920465392256, 1.901297629182468, 2.6120269328999304};
    std::cout << "Predicted Ratios" << '\t' << "Measured Ratios" << '\t' << "(Rpred - Robs) / Robs" << std::endl;
    // Tensor glueball ratio
    std::cout << "TENSOR GLUEBALL SECTOR" << std::endl;
    std::cout << RTGPred[0] << '\t' << RTG_rho[0] << '\t' << (RTGPred[0] - RTG_rho[0]) / RTG_rho[0] << std::endl;
    sum_errors += fabs((RTGPred[0] - RTG_rho[0]) / RTG_rho[0]);
    std::cout << std::endl;
    // Vector meson ratios
    std::cout << "VECTOR MESON NONSINGLET SECTOR" << std::endl;
    for(int i = 0; i < Rrho_rho.size(); i++)
    {
        std::cout << RrhoPred[i] << '\t' << Rrho_rho[i] << '\t' << (RrhoPred[i] - Rrho_rho[i]) / Rrho_rho[i] << std::endl;
        sum_errors += fabs((RrhoPred[i] - Rrho_rho[i]) / Rrho_rho[i]);
    }
    std::cout << std::endl;
    // Axial vector meson ratios
    std::cout << "AXIAL VECTOR MESON NONSINGLET SECTOR" << std::endl;
    for(int i = 0; i < Ra1_rho.size(); i++)
    {
        std::cout << Ra1Pred[i] << '\t' << Ra1_rho[i] << '\t' << (Ra1Pred[i] - Ra1_rho[i]) / Ra1_rho[i] << std::endl;
        sum_errors += fabs((Ra1Pred[i] - Ra1_rho[i]) / Ra1_rho[i]);
    }
    std::cout << std::endl;
    std::cout << "PSEUDOSCALAR MESON NONSINGLET SECTOR" << std::endl;
    for(int i = 0; i < Rpi_rho.size(); i++)
    {
        std::cout << RpiPred[i] << '\t' << Rpi_rho[i] << '\t' << (RpiPred[i] - Rpi_rho[i]) / Rpi_rho[i] << std::endl;
        sum_errors += fabs((RpiPred[i] - Rpi_rho[i]) / Rpi_rho[i]);
    }
    std::cout << std::endl;
    std::cout << "SCALAR MESON NONSINGLET SECTOR" << std::endl;
    for(int i = 0; i < new_Ra0_rho.size(); i++)
    {
        std::cout << Ra0Pred[i] << '\t' << new_Ra0_rho[i] << '\t' << (Ra0Pred[i] - new_Ra0_rho[i]) / new_Ra0_rho[i] << std::endl;
        sum_errors += fabs((Ra0Pred[i] - new_Ra0_rho[i]) / new_Ra0_rho[i]);
    }
    std::cout << std::endl;
    // Singlet Vector meson ratios
    std::cout << "VECTOR MESON SINGLET SECTOR" << std::endl;
    for(int i = 0; i < Romega_rho.size(); i++)
    {
        std::cout << RomegaPred[i] << '\t' << Romega_rho[i] << '\t' << (RomegaPred[i] - Romega_rho[i]) / Romega_rho[i] << std::endl;
        sum_errors += fabs((RomegaPred[i] - Romega_rho[i]) / Romega_rho[i]);
    }
    std::cout << std::endl;
    // Singlet Axial vector meson ratios
    std::cout << "AXIAL VECTOR MESON SINGLET SECTOR" << std::endl;
    for(int i = 0; i < Rf1_rho.size(); i++)
    {
        std::cout << Rf1Pred[i] << '\t' << Rf1_rho[i] << '\t' << (Rf1Pred[i] - Rf1_rho[i]) / Rf1_rho[i] << std::endl;
        sum_errors += fabs((Rf1Pred[i] - Rf1_rho[i]) / Rf1_rho[i]);
    }
    std::cout << std::endl;
    std::cout << "Sum of the relative errors: " << sum_errors << std::endl;

    return 0;
}