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

    invls = 0; fg = 0; fm = 0;
    ag = -7.2104, bg = -7.21643, cg = -6.87905, dg = -11.0486, eg = -15.3461;
    am = -10.8103, bm = -5.35897, cm = -10.8907, dm = -5.58653, em = -9.74851;
    
    // Compute Chebyschev matrices
    chebSetN(400);

    // Get relevant quantities from HVQCD
    //std::cout << hvqcd().QuarkMass() << std::endl;
    std::cout << "umax: " << (hvqcd().u())[0] << std::endl;
    // Setup gluon kernel and compute the Reggeons for t = 0
    GluonKernel gluon(2, {0, ag, bg, cg, dg, eg, 0});
    gluon.computeReggeTrajectories();
    std::vector<Reggeon> reggeons = computeReggeons(gluon, 0.0, 2);

    Reggeon reg = reggeons[0];
    //std::vector<double> z = reg.getWf()[0];
    //std::vector<double> psi = reg.getWf()[1];

    std::cout << "Hard pomeron intercept: " << reg.getJ() << std::endl;

    //std::ofstream myfile;
    //myfile.open("plots/Wavefunction/gluon_wf_0.txt");
    //myfile << "z\tpsi" << std::endl;
    //for(int i = 0; i < z.size(); i++) myfile << z[i] << '\t' << psi[i] << std::endl;
    //myfile.close();


    reg = reggeons[1];
    //z = reg.getWf()[0];
    //psi = reg.getWf()[1];

    std::cout << "Soft pomeron intercept: " << reg.getJ() << std::endl;

    //myfile.open("plots/Wavefunction/gluon_wf_1.txt");
    //myfile << "z\tpsi" << std::endl;
    //for(int i = 0; i < z.size(); i++) myfile << z[i] << '\t' << psi[i] << std::endl;
    //myfile.close();

    // Setup meson kernel and compute the Reggeon for t = 0
    MesonKernel meson(1, {0, am, bm, cm, dm, em, 0});
    meson.computeReggeTrajectories();
    reggeons = computeReggeons(meson, 0.0, 1);

    reg = reggeons[0];
    //z = reg.getWf()[0];
    //psi = reg.getWf()[1];

    std::cout << "Meson intercept: " << reg.getJ() << std::endl;

    //myfile.open("plots/Wavefunction/meson_wf_0.txt");
    //myfile << "z\tpsi" << std::endl;
    //for(int i = 0; i < z.size(); i++) myfile << z[i] << '\t' << psi[i] << std::endl;
    //myfile.close();

    return 0;
}