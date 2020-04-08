#include "schrodinger/common.h"

Point::Point(): x(0), y(0) {}

Point::Point(double xx, double yy): x(xx), y(yy){}


Range::Range(double m0, double m1): eMin(m0), eMax(m1) {}

Mode::Mode(): energy(0), index(-1), wavefunction({}) {}

Mode::Mode(const double e, const std::vector<Point> &f): energy(e), wavefunction(f) {}

Mode::Mode(const double e, const std::vector<Point> &f, const int n): energy(e), index(n), wavefunction(f) {}

void Spectrum::addMode(const Mode &m)
{
    modes.push_back(m);
}

void Spectrum::clear()
{
    modes.clear();
    potential.clear();
}

std::vector<double> Spectrum::getEnergies()
{
    std::vector<double> energies;
    for(int i = 0; i < modes.size(); i++) energies.push_back(modes[i].energy);
    return energies;
}

std::vector<std::vector<Point> > Spectrum::getWavefunctions()
{
    std::vector<std::vector<Point> > wfs;
    for(int i = 0; i < modes.size(); i++) wfs.push_back(modes[i].wavefunction);
    return wfs;
}