#ifndef F2PHOTON_H
#define F2PHOTON_H

#include <vector>
#include "DeepInelasticScattering.h"

class F2Photon: public DeepInelasticScattering
{
    public :
        F2Photon(std::string data_path = "expdata/F2Photon/F2Photon_data.txt");
        double IzN(const std::vector<double> &kin, const Reggeon &reg);
};

#endif