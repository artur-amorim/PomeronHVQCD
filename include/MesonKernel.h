#ifndef MESON_KERNEL_H
#define MESON_KERNEL_H

#include <vector>
#include "Kernel.h"

class MesonKernel : public Kernel{
    public:
        // class constructor
        MesonKernel(const int nreg, const std::vector<double> &pars);
        // Computes the potential of the gluon kernel
        std::vector<double> computePotential(const double J) const;
};
#endif