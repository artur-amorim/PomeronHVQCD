#ifndef KERNEL_H
#define KERNEL_H

#include <vector>
#include <string>
#include "Reggeon.h"
#include "methods/interpolation/Spline_Interp.hpp"

class Kernel
{
    private :
        std::string name;                                         // Name of the Kernel
        int nTrajectories;                                        // Number of desired Regge trajectories
        std::vector<Spline_Interp<double> > trajectories;               // Vector containing the reggeons
        std::vector<double> kernelPars;                                 // Parameters that characterize the potential
    protected:
        void copy(const Kernel &rhs);
    public :
        // Kernel constructors
        Kernel();                                                  
        Kernel(const std::string nname, const int ntraj, const std::vector<double> &kerpars);
        Kernel(const Kernel &kernel);                                                            // Copy constructor
        std::string Name() const;                                                                // Getter of name
        int NTrajectories() const ;                                                              // Getter of nTrajectories
        Spline_Interp<double> Trajectory(const int n) const;                                     // Returns the nth trajectory
        std::vector<double> KernelPars() const;                                                  // Getter of kernelPars
        void computeReggeTrajectories(const std::vector<double> &pars = {});                     // Computes the Regge Trajectories
        void setKernelPars(const std::vector<double> &pars);                                     // Setter of kernelPars
        virtual std::vector<double> computePotential(const double J) const;                      // Compute the potential vector of a kernel
        Kernel& operator= (const Kernel &rhs);                                                   // Assignment operator
        ~Kernel();                                                                               // Kernel destructor
};

std::vector<Reggeon> computeReggeons(const Kernel &kernel, const double t, const int n);          // Function that computes the reggeons.

#endif