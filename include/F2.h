#ifndef F2_H
#define F2_H

#include <vector>
#include "DeepInelasticScattering.h"

class F2 : public DeepInelasticScattering {
    public :
        F2(std::string file_path = "expdata/F2_data.txt");                                         // Constructor
        double IzN(const std::vector<double> &kin, const Reggeon &reg);
};

#endif