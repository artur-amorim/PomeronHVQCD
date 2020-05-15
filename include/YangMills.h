#ifndef YANGMILLS_H
#define YANGMILLS_H

#include <vector>
#include "background.h"

class YangMills: public Background
{
    private:
        std::vector<double> d3Phis;
        // Computes d3Phi/dA3 in Yang-Mills theory
        double d3PhiYM(const double q, const double phi);
        // Yang Mills EOMs
        void eom(const state &X, state &dXdA, const double A);
        // defines how to save the computed background fields for Yang-Mills
        void observer(const state &X, const double A);
        // prepare background for further computations
        void finalizeBackground(const double AIR, const double AUV);
        void solveRaw(const double AIR, const double AUV);
    public:
        // YangMills constructor
        YangMills(const double ssc = 3.0, const double VVgIR = 2.05);
        // YangMills copy constructor
        YangMills(const YangMills &ym);
        // YangMills assignment operator
        YangMills& operator=(const YangMills &rhs); 
        // YangMills destructor
        ~YangMills();
        // Get d3Phi values
        std::vector<double> d3Phi() const;
};

// Computes the scalar 0^{++} glueball potential
std::vector<double> computeV0GPotential(const YangMills &ym);


// Computes the Yang-Mills theory spectrum
void computeYangMillsSpectrum(const YangMills &ym, const int n_tensor = 1, const int n_scalar = 2);

// Computes the Yang-Mills theory glueball ratios
void computeYangMillsRatios(const YangMills &ym);

const double R20 = 1.46;
const double R00 = 1.87;

// Fits the parameters to the Yang-Mills ratios
void fitYangMills(const double sc_guess, const double VgIR_guess);


#endif