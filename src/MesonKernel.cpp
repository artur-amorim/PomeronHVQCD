#include <cmath>
#include "HolographicVQCD.h"
#include "MesonKernel.h"
#include "methods/vectorOperators.hpp"

// GluonKernel constructor
MesonKernel::MesonKernel(const int nreg, const std::vector<double> &pars):
    Kernel("meson", nreg, pars) {}

std::vector<double> MesonKernel::computePotential(const double J) const
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
    // model = IHQCD
    std::vector<double> VSch = hvqcd().U1() + (J*J-1) * hvqcd().e2A();
    VSch = VSch + (J-1) * ( (1.0 + invls * invls) * pow(hvqcd().dAstring(), 2) - hvqcd().d2Astring() / pow(hvqcd().G(), 2) );
    //VSch = VSch + (J-1) * (a * hvqcd().aF() + b * hvqcd().bF() + c * hvqcd().cF() + d * hvqcd().dF() + e * hvqcd().eF());
    // VSch's values are ordered from IR to UV
    // We want from UV to IR, so we reverse them
    std::reverse(std::begin(VSch), std::end(VSch));
    return VSch;
}