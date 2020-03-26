#ifndef GLUON_KERNEL_H
#define GLUON_KERNEL_H

#include <vector>
#include "Kernel.h"

class GluonKernel : public Kernel{
    public:
        // class constructor
        GluonKernel(const int nreg, const std::vector<double> &pars);
        // Computes the potential of the gluon kernel
        std::vector<double> computePotential(const double J) const;
};
#endif