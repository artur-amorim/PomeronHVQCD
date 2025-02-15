#ifndef SCHRODINGER_H
#define SCHRODINGER_H

#include <vector>

class SolvSpec;

struct List
{
  std::vector<double> Es;
  std::vector< std::vector< std::vector<double> > > wfs;
  List(const std::vector<double> &es = {}, const std::vector< std::vector< std::vector<double> > > &Wfs = {})
  {
    Es = es;
    wfs = Wfs;
  }
};

SolvSpec* setSchroMethod(std::string method);
List computeSpectrum(const std::vector<double> &px , const std::vector<double> &py,
                     int nEigen = 3, std::string method = "cheb",
                     double dE = 0.1, double tol = 1e-9);
;
#endif