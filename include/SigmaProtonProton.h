#ifndef SIGPP_H
#define SIGPP_H

#include <vector>
#include "Sigma.h"

class SigmaProtonProton: public Sigma
{
    public :
        SigmaProtonProton(std::string data_path = "");
        double IzN(const std::vector<double> &kin, const Reggeon &reg);
};

#endif