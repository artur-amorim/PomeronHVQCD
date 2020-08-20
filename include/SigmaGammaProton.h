#ifndef SIGGAMMAP_H
#define SIGGAMMAP_H

#include <vector>
#include "Sigma.h"

class SigmaGammaProton: public Sigma
{
    public :
        SigmaGammaProton(std::string data_path = "");
};

#endif