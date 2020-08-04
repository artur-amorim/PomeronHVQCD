#include <iostream>
#include <fstream>
#include "HolographicVQCD.h"
#include "methods/interpolation/Poly_Interp.hpp"

int main()
{
    std::vector<double> u = hvqcd().u();
    std::vector<double> Vfw2fac(u.size()), MesonPotFac(u.size());
    // e^(-7/3 \Phi) V_f w_s^2 = e^(\Phi / 3) V_f w^2 in terms of the Einstein frame background potentials
    // e^(-10/3 \Phi) V_f w_s^2 = e^(-2 \Phi / 3) V_f w^2 in terms of the Einstein frame background potentials
    std::vector<double> Phis = hvqcd().Phi(), taus = hvqcd().tau();
    for(int i = 0; i < u.size(); i++) 
    {
        Vfw2fac[i] = std::exp(Phis[i] / 3) * hvqcd().Vf(Phis[i], taus[i]) * std::pow(hvqcd().w(Phis[i]),2);
        MesonPotFac[i] = std::sqrt(std::exp(-2*Phis[i]/3) * hvqcd().Vf(Phis[i], taus[i]) * std::pow(hvqcd().w(Phis[i]),2));
    }
    // Now we reverse u, Vfw2fac and MesonPotFac because we want them from the UV to the IR
    std::reverse(std::begin(u), std::end(u));
    std::reverse(std::begin(Vfw2fac), std::end(Vfw2fac));
    std::reverse(std::begin(MesonPotFac), std::end(MesonPotFac));
    Poly_Interp<double> potFactor(u, Vfw2fac, 4), MesonPotFactor(u, MesonPotFac, 4);
    
    std::ofstream file;
    file.open("PotProfiles.txt");
    file << "u\tGluePotFac\tMesonPotFac" << std::endl;
    for(int i = 0; i < u.size(); i++) file << u[i] << '\t' << Vfw2fac[i] << '\t' << MesonPotFac[i] << std::endl;
    file.close();

    file.open("PotInterpolationProfiles.txt");
    file << "u\tGluePotFac\tMesonPotFac" << std::endl;
    for(double u = 0.001; u <= 10; u += 0.01) file << u << '\t' << potFactor.interp(u) << '\t' << MesonPotFactor.interp(u) << std::endl;
    file.close();
    
    return 0;
}