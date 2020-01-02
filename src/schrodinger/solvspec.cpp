#include <iostream>
#include <fstream>

#include "schrodinger/solvspec.h"

SolvSpec::SolvSpec(){}

void SolvSpec::setPotential(const std::vector<double> &XX, const std::vector<double> &VV){
	if (XX.size() != VV.size())
	{
		std::cout << "X and V must have the same size" << std::endl;
		throw "error";
	}
	if(XX.size() > 1) 
	{
    	X = XX;
		PotVals = VV;
  	} 
	else 
	{
		std::cout << "Give more points to define the potential." << std::endl;
		throw "error"; 
  }

  nPoints = XX.size() ;
  xMin    = XX.front() ;
  xMax    = XX.back() ;
  h = (xMax - xMin) / nPoints;
  // Now we define the Poly_Interp object that will be used to compute the potential at any x
  potFunc = Poly_Interp<double>(XX, VV, 4);
}

std::vector<Point> SolvSpec::getPotential()
{
  // Returns the potential attribute
	std::vector<Point> v;
	for(int i = 0; i < nPoints; i++)
	{
		v.push_back(Point(X[i], PotVals[i]));
	}
	return v;
}

void SolvSpec::findSpectrum(int nEigen){}

Spectrum SolvSpec::getSpectrum()
{
	// Returns the spectrum attribute
	return spectrum;
}

void SolvSpec::savePotential() 
{
  	std::ofstream f("potential.dat");
	for (double x = xMin; x < xMax; x += 0.1) 
	{
		f << x << " " << potFunc.interp(x) << std::endl;
	}
	f.close();
}

SolvSpec::~SolvSpec(){}