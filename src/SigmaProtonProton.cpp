#include <iostream>
#include <cmath>
#include "SigmaProtonProton.h"

// Convert sigma(p p -> X) values to the appropriate units.
const double mb_to_GEVMINUS2 = 10 / 3.894 ;
const double mrho_U_mrho_GeV = 4.3669;
const double mb_to_UMINUS2 = mb_to_GEVMINUS2 / std::pow(mrho_U_mrho_GeV,2);

SigmaProtonProton::SigmaProtonProton(std::string data_path): Sigma(data_path, mb_to_GEVMINUS2) 
{
    std::cout << "Loaded sigma(p p -> X) data." << std::endl;
}

double SigmaProtonProton::IzN(const std::vector<double> &kin, const Reggeon &reg)
{
    return 1;
}