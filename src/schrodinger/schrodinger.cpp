#include <string>
#include "schrodinger/schrodinger.h"


SolvSpec* setSchroMethod(std::string method) 
{
  if(method == "numerov") {
    return new Numerov();
  }
  // use chebyshev as default method
  return new ChebSpec();
}

std::vector< std::vector<double> > getPotential(SolvSpec* n) 
{
  std::vector<Point> p = n->getPotential();
  int length = p.size();

  std::vector<double> x, y;
  for(int i = 0; i < length; i++) 
  {
    x.push_back(p[i].x);
    y.push_back(p[i].y);
  }

  std::vector< std::vector<double> > pot{x, y};
  return pot;
}

std::vector<double> getEnergies(SolvSpec* n) 
{
  std::vector<double> energies = n->getSpectrum().getEnergies();
  return energies;
}

std::vector< std::vector< std::vector<double> > > getWavefunctions(SolvSpec* n) 
{
  std::vector< std::vector< std::vector<double> > > wfs;
  std::vector<std::vector<Point> > WFs = n->getSpectrum().getWavefunctions();

  for(int i = 0; i < WFs.size(); i++) 
  {
    int length = WFs[i].size();
    std::vector<double> x, y;
    for(int j = 0; j < length; j++) 
    {
      x.push_back(WFs[i][j].x);
      y.push_back(WFs[i][j].y);
    }
    std::vector< std::vector<double> > wf{x, y};
    wfs.push_back(wf);
  }
  return wfs;
}

List computeSpectrum(const std::vector<double> &px , const std::vector<double> &py,
                     int nEigen, std::string method,
                     double dE, double tol) 
{
  
  SolvSpec* n = setSchroMethod(method);
  
  if(px.size() != py.size()) 
  {
    std::cout << px.size() << '\t' << py.size() << std::endl;
    std::cout << "Please pass two columns with the same size for the potential" << std::endl;
    throw "error";
  }

  n->setPotential(px, py);
  n->dEmin = dE;
  n->tol = tol;
  n->findSpectrum(nEigen);
  std::vector<double> energies = getEnergies(n) ;
  std::vector<std::vector<std::vector<double> > > wavefuncs = getWavefunctions(n) ;
  // Free memory
  delete n ;
  return List(energies, wavefuncs) ;
}