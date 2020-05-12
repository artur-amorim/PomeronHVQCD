#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <vector>
#include <cmath>
#include <functional>

struct Point
{
    double x, y;
    Point();
    Point(double xx, double yy);
};

struct Range
{
	double eMin, eMax;
	Range(double m0, double m1);
};

struct Mode
{
	double energy;
  	int index;
	std::vector<Point> wavefunction;
	Mode();
	Mode(const double e, const std::vector<Point> &f);
	Mode(const double e, const std::vector<Point> &f, const int n);
};

struct Spectrum
{
	std::vector<Mode>  modes;
	std::vector<std::vector<double> > potential;

	void addMode(const Mode &m);
	void clear();
	std::vector<double> getEnergies();
  	std::vector<std::vector<Point> > getWavefunctions();
};

#endif