#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <vector>
#include <cmath>
#include <functional>

template<typename T>
inline T SIGN(const T& a, const T &b)
{
	return (b) >= 0.0 ? std::fabs(a) : -std::fabs(a) ;
}

struct Point
{
    double x, y;
    Point() 
	{
        x = 0;
        y = 0;
    }
    Point(double xx, double yy)
	{
        x = xx;
        y = yy;
    }
};

struct Range
{
	double eMin, eMax;
	Range(double m0, double m1)
	{
	    eMin = m0;
	    eMax = m1;
	}
};

struct Mode
{
	double energy;
  	int index;
	std::vector<Point> wavefunction;
	Mode() 
	{
		// default constructor for non good modes
    	index = -1;
	}
	Mode(double e, std::vector<Point> f)
	{
		energy = e;
		wavefunction = f;
	}
	Mode(double e, std::vector<Point> f, int n)
	{
		energy = e;
		wavefunction = f;
		index = n;
	}
};

struct Spectrum
{
	std::vector<Mode>  modes;
	std::vector<std::vector<double> > potential;

	void addMode(Mode m)
	{
   		modes.push_back(m);
	}

	void clear() 
	{
   		modes.clear();
	  	potential.clear();
 	}
	
	std::vector<double> getEnergies()
	{
		std::vector<double> energies;
		for(int i = 0; i < modes.size(); i++)
		{
      		energies.push_back(modes[i].energy);
		}
    	return energies;
	}

  std::vector<std::vector<Point> > getWavefunctions()
  {
    std::vector<std::vector<Point> > wfs;
    for(int i = 0; i < modes.size(); i++)
	{
      wfs.push_back(modes[i].wavefunction);
	}
    return wfs;
  }
};

// Van Wijngaarden–Dekker–Brent Method for finding root, from NR
double zbrent(std::function<double(double)>& func, double x1, double x2, double tol, bool silent = false);

double bisection(std::function<double(double)> diffFunc, double min, double max);

#endif