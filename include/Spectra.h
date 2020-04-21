#ifndef SPECTRA_H
#define SPECTRA_H

#include <vector>
#include "Reggeon.h"

class Spectra
{
    private :
        double t;
        std::vector<Reggeon> reggeons;
        void copy(const Spectra &rhs);
    public :
        // Default constructor of Spectra class
        Spectra(const double tval = 0);
        // Spectra class constructor
        Spectra(const double &tval, const std::vector<Reggeon> &regs);
        // Spectra class copy constructor
        Spectra(const Spectra &spec);
        // Setter of ts
        void setT(const double tval);
        // Setter of kernels
        void setReggeons(const std::vector<Reggeon> &regs);
        // Getter of t
        double getT() const;
        // Getter of kernels
        std::vector<Reggeon> getReggeons() const;
        // Assignment operator
        Spectra& operator= (const Spectra &rhs);
        // operator<
        bool operator< (const Spectra &rhs);
        // operator>
        bool operator> (const Spectra &rhs);
        // Destructor
        ~Spectra();
};
#endif