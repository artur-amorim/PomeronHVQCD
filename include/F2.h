#ifndef F2_H
#define F2_H

#include <vector>
#include "DeepInelasticScattering.h"

class F2 : public DeepInelasticScattering {
    public :
        F2(const bool rrsslog = false, std::string file_path = "expdata/F2_data.txt");                                         // Constructor
        double IzN(const std::vector<double> &kin, const Reggeon &reg);
};

class F2IzNIntegrand: public IzNIntegrand
{
    public:
        F2IzNIntegrand(const Poly_Interp<double> &f1, const Poly_Interp<double> &f2, const U1NNMode &f3, const Poly_Interp<double> &f4);
        double operator()(const double x);
};

#endif