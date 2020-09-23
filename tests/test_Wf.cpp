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

    invls = -0.00329973; ag = -4.42448; bg = -11.9668; cg = 2.00837; dg = 2.75436; eg = -0.229314; fg = 4.79135;
    am = 4.5345; bm = -4.46109; cm = 10.2159; dm = 0.444645; em = 6.65785; fm = 1.16794;
    
    // Compute Chebyschev matrices
    chebSetN(400);

    // Get relevant quantities from HVQCD
    //std::cout << hvqcd().QuarkMass() << std::endl;
    std::cout << "umax: " << (hvqcd().u())[0] << std::endl;
    // Setup gluon kernel and compute the Reggeons for t = 0
    GluonKernel gluon(2, {invls, ag, bg, cg, dg, eg, fg});
    gluon.computeReggeTrajectories();
    std::vector<Reggeon> reggeons = computeReggeons(gluon, 0.0, 2);

    Reggeon reg = reggeons[0];
    std::vector<double> z = reg.getWf()[0];
    std::vector<double> psi = reg.getWf()[1];

    std::cout << "Hard pomeron intercept: " << reg.getJ() << std::endl;

    std::ofstream myfile;
    myfile.open("plots/Wavefunction/gluon_wf_0.txt");
    myfile << "z\tpsi" << std::endl;
    for(int i = 0; i < z.size(); i++) myfile << z[i] << '\t' << psi[i] << std::endl;
    myfile.close();


    reg = reggeons[1];
    z = reg.getWf()[0];
    psi = reg.getWf()[1];

    std::cout << "Soft pomeron intercept: " << reg.getJ() << std::endl;

    myfile.open("plots/Wavefunction/gluon_wf_1.txt");
    myfile << "z\tpsi" << std::endl;
    for(int i = 0; i < z.size(); i++) myfile << z[i] << '\t' << psi[i] << std::endl;
    myfile.close();

    // Setup meson kernel and compute the Reggeon for t = 0
    MesonKernel meson(1, {invls, am, bm, cm, dm, em, fm});
    meson.computeReggeTrajectories();
    reggeons = computeReggeons(meson, 0.0, 1);

    reg = reggeons[0];
    z = reg.getWf()[0];
    psi = reg.getWf()[1];

    std::cout << "Meson intercept: " << reg.getJ() << std::endl;

    myfile.open("plots/Wavefunction/meson_wf_0.txt");
    myfile << "z\tpsi" << std::endl;
    for(int i = 0; i < z.size(); i++) myfile << z[i] << '\t' << psi[i] << std::endl;
    myfile.close();

    return 0;
}