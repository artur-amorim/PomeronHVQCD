#ifndef SOLVSPEC_H
#define SOLVSPEC_H

#include <vector>
#include "common.h"
#include "../methods/interpolation/Poly_Interp.hpp"

class SolvSpec {
  protected:
  	std::vector<double> X, PotVals;
  	Spectrum spectrum;
  public:
  	double xMin, xMax, tol;
    // Attributes relevant for the Numerov Method
  	double h, dEmin;
    int nPoints;
    // Object with interp method that computes the potential at arbitrary x
    Poly_Interp<double> potFunc;
    // Class constructor
  	SolvSpec();
    // Setter of potential
  	virtual void setPotential(const std::vector<double> &XX, const std::vector<double> &VV);
    // Getter of potential
  	std::vector<Point> getPotential();
    // Computes the spectrum using a given method
  	virtual void findSpectrum(int nEigen);
    // Get the spectrum
  	Spectrum getSpectrum();
    // Save the potential in a file
  	void savePotential();
    // Class destructor
	virtual ~SolvSpec();
};

#endif