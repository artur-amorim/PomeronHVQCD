#include <string>
#include <memory>
#include "schrodinger/schrodinger.h"
#include "schrodinger/chebspec.h"
#include "schrodinger/numerov.h"
#include "schrodinger/matrixNumerov.h"


SolvSpec* setSchroMethod(std::string method) 
{
  if(method == "numerov") return new Numerov();
  else if (method == "matrix_numerov") return new MatrixNumerov();
  // use chebyshev as default method
  return new ChebSpec();
}

List computeSpectrum(const std::vector<double> &px , const std::vector<double> &py,
                     int nEigen, std::string method,
                     double dE, double tol) 
{
  
  // Create smart pointer of SolvSpec. No need to worry about deleting the pointer n;
  std::unique_ptr<SolvSpec> n(setSchroMethod(method));
  
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
  // Get the energies
  std::vector<double> energies = n->getSpectrum().getEnergies();
  // Get the wavefunctions
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
  return List(energies, wfs) ;
}