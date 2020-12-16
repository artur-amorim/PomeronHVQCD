#include "HolographicVQCD.h"
#include "GluonKernel.h"
#include "methods/vectorOperators.hpp"

// GluonKernel constructor
GluonKernel::GluonKernel(const int nreg, const std::vector<double> &pars):
    Kernel("gluon", nreg, pars) {}

std::vector<double> GluonKernel::computePotential(const double J) const
{
    // Computes the potential values given J
    // Get the kernel parameters kernelPars = {invls, a, b ,c , d, e, f}
    std::vector<double> pars = this->KernelPars();
    const double invls = pars[0];
    const double a = pars[1];
    const double b = pars[2];
    const double c = pars[3];
    const double d = pars[4];
    const double e = pars[5];
    const double f = pars[6];
    // model = HVQCD
    std::vector<double> VSch = hvqcd().U2() + (J-2) * 2.0 * invls * invls * hvqcd().e2Astring() ;
    //std::vector<double> VSch = hvqcd().U2() + (J-2) * ( 2.0 * invls * invls * hvqcd().e2Astring() * (1.0 + f / hvqcd().l1_2()) + (J+2) * hvqcd().e2A()  + a ///* hvqcd().aF() + b * hvqcd().bF()  + c * hvqcd().cF() + d * hvqcd().dF() + e * hvqcd().eF()) ;
    // VSch's values are ordered from IR to UV
    // We want from UV to IR, so we reverse them
    std::reverse(std::begin(VSch), std::end(VSch));
    return VSch;
}