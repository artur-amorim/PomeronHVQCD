#include <iostream>
#include "SigmaGammaGamma.h"

const double nb_to_GEVMINUS2 = 1.0 / (3.894e5) ;
const double mrho_U_mrho_GeV = 4.3669;
const double nb_to_UMINUS2 = nb_to_GEVMINUS2 / std::pow(mrho_U_mrho_GeV,2);

SigmaGammaGamma::SigmaGammaGamma(std::string data_path): Sigma(data_path, nb_to_UMINUS2)
{
    std::cout << "Loaded sigma(gamma gamma -> X) data." << std::endl;
}