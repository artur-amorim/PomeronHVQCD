#include <iostream>
#include <fstream>
#include <vector>
#include "HolographicVQCD.h"
#include "Reggeon.h"
#include "GluonKernel.h"
#include "schrodinger/schrodinger.h"
#include "schrodinger/chebspec.h"

int main()
{
    // Compute Chebyschev matrices
    chebSetN(400);

    // Get relevant quantities from HVQCD
    std::cout << hvqcd().QuarkMass() << std::endl;
    std::vector<double> zvals = hvqcd().z();
    std::reverse(std::begin(zvals), std::end(zvals));
    
    // Setup gluon kernel and compute the Reggeons for t = 0
    std::vector<double> gluon_pars = {0.232303, 7.33966, 0.696812, 1.15569, 0.569428, 0.940377, -0.27283};
    GluonKernel gluon(4, gluon_pars);

    gluon.computeReggeTrajectories();
    std::vector<Reggeon> reggeons = computeReggeons(gluon, 0.0, 4);

    Reggeon reg = reggeons[0];
    std::vector<double> z = reg.getWf()[0];
    std::vector<double> psi = reg.getWf()[1];
    std::vector<double> potential = gluon.computePotential(reg.getJ());

    std::ofstream myfile;
    myfile.open("plots/Wavefunction/wf_0.txt");
    myfile << "z\tpsi" << std::endl;
    for(int i = 0; i < z.size(); i++) myfile << z[i] << '\t' << psi[i] << std::endl;
    myfile.close();

    myfile.open("plots/Wavefunction/potential_0.txt");
    myfile << "z\tpotential" << std::endl;
    for(int i = 0; i < zvals.size(); i++) myfile << zvals[i] << '\t' << potential[i] << std::endl;
    myfile.close();

    reg = reggeons[1];
    z = reg.getWf()[0];
    psi = reg.getWf()[1];
    potential = gluon.computePotential(reg.getJ());

    myfile.open("plots/Wavefunction/wf_1.txt");
    myfile << "z\tpsi" << std::endl;
    for(int i = 0; i < z.size(); i++) myfile << z[i] << '\t' << psi[i] << std::endl;
    myfile.close();

    myfile.open("plots/Wavefunction/potential_1.txt");
    myfile << "z\tpotential" << std::endl;
    for(int i = 0; i < zvals.size(); i++) myfile << zvals[i] << '\t' << potential[i] << std::endl;
    myfile.close();

    reg = reggeons[2];
    z = reg.getWf()[0];
    psi = reg.getWf()[1];
    potential = gluon.computePotential(reg.getJ());

    myfile.open("plots/Wavefunction/wf_2.txt");
    myfile << "z\tpsi" << std::endl;
    for(int i = 0; i < z.size(); i++) myfile << z[i] << '\t' << psi[i] << std::endl;
    myfile.close();

    myfile.open("plots/Wavefunction/potential_2.txt");
    myfile << "z\tpotential" << std::endl;
    for(int i = 0; i < zvals.size(); i++) myfile << zvals[i] << '\t' << potential[i] << std::endl;
    myfile.close();

    reg = reggeons[3];
    z = reg.getWf()[0];
    psi = reg.getWf()[1];
    potential = gluon.computePotential(reg.getJ());

    myfile.open("plots/Wavefunction/wf_3.txt");
    myfile << "z\tpsi" << std::endl;
    for(int i = 0; i < z.size(); i++) myfile << z[i] << '\t' << psi[i] << std::endl;
    myfile.close();

    myfile.open("plots/Wavefunction/potential_3.txt");
    myfile << "z\tpotential" << std::endl;
    for(int i = 0; i < zvals.size(); i++) myfile << zvals[i] << '\t' << potential[i] << std::endl;
    myfile.close();

    potential = gluon.computePotential(2);
    std::vector<double> glueball_potential = computeV2GPotential(hvqcd());
    std::reverse(std::begin(glueball_potential), std::end(glueball_potential));

    myfile.open("plots/Wavefunction/gluon_pot_J_2_vs_glueball_potential.txt");
    myfile << "z\tglue\tglueball" << std::endl;
    for(int i = 0; i < zvals.size(); i++) myfile << zvals[i] << '\t' << potential[i] << '\t' << glueball_potential[i] << std::endl;
    myfile.close();

    std::vector<double> aF = hvqcd().aF();
    std::reverse(std::begin(aF), std::end(aF));
    myfile.open("plots/Wavefunction/aF_vs_z.txt");
    myfile << "z\taF" << std::endl;
    for(int i = 0; i < zvals.size(); i++) myfile << zvals[i] << '\t' << aF[i] << std::endl;
    myfile.close();

    std::vector<double> bF = hvqcd().bF();
    std::reverse(std::begin(bF), std::end(bF));
    myfile.open("plots/Wavefunction/bF_vs_z.txt");
    myfile << "z\tbF" << std::endl;
    for(int i = 0; i < zvals.size(); i++) myfile << zvals[i] << '\t' << bF[i] << std::endl;
    myfile.close();

    std::vector<double> cF = hvqcd().aF();
    std::reverse(std::begin(cF), std::end(cF));
    myfile.open("plots/Wavefunction/cF_vs_z.txt");
    myfile << "z\tcF" << std::endl;
    for(int i = 0; i < zvals.size(); i++) myfile << zvals[i] << '\t' << cF[i] << std::endl;
    myfile.close();

    std::vector<double> dF = hvqcd().dF();
    std::reverse(std::begin(dF), std::end(dF));
    myfile.open("plots/Wavefunction/dF_vs_z.txt");
    myfile << "z\tdF" << std::endl;
    for(int i = 0; i < zvals.size(); i++) myfile << zvals[i] << '\t' << dF[i] << std::endl;
    myfile.close();

    std::vector<double> eF = hvqcd().eF();
    std::reverse(std::begin(eF), std::end(eF));
    myfile.open("plots/Wavefunction/eF_vs_z.txt");
    myfile << "z\teF" << std::endl;
    for(int i = 0; i < zvals.size(); i++) myfile << zvals[i] << '\t' << eF[i] << std::endl;
    myfile.close();

    std::vector<double> e2A = hvqcd().e2A();
    std::reverse(std::begin(e2A), std::end(e2A));
    myfile.open("plots/Wavefunction/e2A_vs_z.txt");
    myfile << "z\te2A" << std::endl;
    for(int i = 0; i < zvals.size(); i++) myfile << zvals[i] << '\t' << e2A[i] << std::endl;
    myfile.close();

    std::vector<double> e2Astring = hvqcd().e2Astring();
    std::reverse(std::begin(e2Astring), std::end(e2Astring));
    myfile.open("plots/Wavefunction/e2Astring_vs_z.txt");
    myfile << "z\te2Astring" << std::endl;
    for(int i = 0; i < zvals.size(); i++) myfile << zvals[i] << '\t' << e2Astring[i] << std::endl;
    myfile.close();

    std::vector<double> l1_2 = hvqcd().l1_2();
    std::reverse(std::begin(l1_2), std::end(l1_2));
    myfile.open("plots/Wavefunction/l1_2_vs_z.txt");
    myfile << "z\tl1_2" << std::endl;
    for(int i = 0; i < zvals.size(); i++) myfile << zvals[i] << '\t' << l1_2[i] << std::endl;
    myfile.close();




    return 0;
}