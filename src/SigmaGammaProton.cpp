#include <iostream>
#include "SigmaGammaProton.h"

const double mub_to_GEVMINUS2 = 1.0 / (3.894e2) ;
const double mrho_U_mrho_GeV = 4.3669;
const double mub_to_UMINUS2 = mub_to_GEVMINUS2 / std::pow(mrho_U_mrho_GeV,2);


SigmaGammaProton::SigmaGammaProton(std::string data_path): Sigma(data_path, mub_to_UMINUS2) 
{
    std::cout << "Loaded sigma(gamma p -> X) data." << std::endl;
}