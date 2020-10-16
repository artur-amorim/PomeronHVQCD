#include <iostream>
#include <fstream>
#include <vector>
#include "HolographicVQCD.h"
#include "Reggeon.h"
#include "GluonKernel.h"
#include "MesonKernel.h"
#include "schrodinger/schrodinger.h"
#include "schrodinger/chebspec.h"

int main()
{
    double invls, ag, bg, cg, dg, eg, fg;
    double am, bm, cm, dm, em, fm;

    invls = 0; ag = -0.952001; bg = 12.8395; cg = -3.62051; dg = 4.85515; eg = 1.0048; fg = 0;
    am = -3.74852; bm = 4.81612; cm = -9.96161; dm = 0.247632; em = 0.115675; fm = 0;
    
    // Compute Chebyschev matrices
    chebSetN(400);

    // Get relevant quantities from HVQCD
    //std::cout << hvqcd().QuarkMass() << std::endl;

    // Create interpolation function between u and z
    std::vector<double> z = hvqcd().z(), u = hvqcd().u();
    std::cout << "umax: " << u[0] << std::endl;
    // Reverse because we want function from UV to IRß
    std::reverse(std::begin(z), std::end(z));
    std::reverse(std::begin(u), std::end(u));
    Poly_Interp<double> ufunc(z, u, 4);

    // Setup gluon kernel and compute the Reggeons for t = 0
    GluonKernel gluon(4, {0, ag, bg, cg, dg, eg, 0});
    gluon.computeReggeTrajectories();
    std::vector<Reggeon> reggeons = computeReggeons(gluon, 0.0, 4);

    Reggeon reg = reggeons[0];
    z = reg.getWf()[0];
    for(int i = 0; i < z.size(); i++) z[i] = ufunc.interp(z[i]);
    std::vector<double> psi = reg.getWf()[1];

    std::cout << "Hard pomeron intercept: " << reg.getJ() << std::endl;

    std::ofstream myfile;
    myfile.open("plots/Wavefunction/gluon_wf_0.txt");
    myfile << "u\tpsi" << std::endl;
    for(int i = 0; i < z.size(); i++) myfile << z[i] << '\t' << psi[i] << std::endl;
    myfile.close();


    reg = reggeons[1];
    z = reg.getWf()[0];
    for(int i = 0; i < z.size(); i++) z[i] = ufunc.interp(z[i]);
    psi = reg.getWf()[1];

    std::cout << "Soft pomeron intercept: " << reg.getJ() << std::endl;

    myfile.open("plots/Wavefunction/gluon_wf_1.txt");
    myfile << "u\tpsi" << std::endl;
    for(int i = 0; i < z.size(); i++) myfile << z[i] << '\t' << psi[i] << std::endl;
    myfile.close();

    reg = reggeons[2];
    z = reg.getWf()[0];
    for(int i = 0; i < z.size(); i++) z[i] = ufunc.interp(z[i]);
    psi = reg.getWf()[1];

    std::cout << "3rd intercept: " << reg.getJ() << std::endl;

    myfile.open("plots/Wavefunction/gluon_wf_2.txt");
    myfile << "u\tpsi" << std::endl;
    for(int i = 0; i < z.size(); i++) myfile << z[i] << '\t' << psi[i] << std::endl;
    myfile.close();


    reg = reggeons[3];
    z = reg.getWf()[0];
    for(int i = 0; i < z.size(); i++) z[i] = ufunc.interp(z[i]);
    psi = reg.getWf()[1];

    std::cout << "4th intercept: " << reg.getJ() << std::endl;

    myfile.open("plots/Wavefunction/gluon_wf_3.txt");
    myfile << "u\tpsi" << std::endl;
    for(int i = 0; i < z.size(); i++) myfile << z[i] << '\t' << psi[i] << std::endl;
    myfile.close();


    // Setup meson kernel and compute the Reggeon for t = 0
    MesonKernel meson(2, {0, am, bm, cm, dm, em, 0});
    meson.computeReggeTrajectories();
    reggeons = computeReggeons(meson, 0.0, 2);

    reg = reggeons[0];
    z = reg.getWf()[0];
    psi = reg.getWf()[1];

    std::cout << "Meson intercept: " << reg.getJ() << std::endl;

    myfile.open("plots/Wavefunction/meson_wf_0.txt");
    myfile << "u\tpsi" << std::endl;
    for(int i = 0; i < z.size(); i++) myfile << z[i] << '\t' << psi[i] << std::endl;
    myfile.close();
    
    reg = reggeons[1];
    z = reg.getWf()[0];
    psi = reg.getWf()[1];

    std::cout << "Meson intercept: " << reg.getJ() << std::endl;

    myfile.open("plots/Wavefunction/meson_wf_1.txt");
    myfile << "u\tpsi" << std::endl;
    for(int i = 0; i < z.size(); i++) myfile << z[i] << '\t' << psi[i] << std::endl;
    myfile.close();

    return 0;
}