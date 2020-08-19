#ifndef FL_H
#define FL_H

#include <vector>
#include "DeepInelasticScattering.h"

class FL : public DeepInelasticScattering
{
    public :
        FL(std::string file_path = "expdata/FL_data.txt");                                           // Constructor
        double IzN(const std::vector<double> &kin, const Reggeon &reg);
};

class FLIzNIntegrand: public IzNIntegrand
{
    public:
        FLIzNIntegrand(const Poly_Interp<double> &f1, const Poly_Interp<double> &f2, const U1NNMode &f3, const Poly_Interp<double> &f4);
        double operator()(const double x);
};

#endif