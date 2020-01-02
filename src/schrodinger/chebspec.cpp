#include <algorithm>
#include "schrodinger/chebspec.h"

// Definition of the struct EigenPair
EigenPair::EigenPair(const double eigval, const arma::vec & eigvec): eigen_value(eigval), eigen_vector(eigvec){}

bool compEigenPair(const EigenPair &lhs, const EigenPair &rhs)
{
  return lhs.eigen_value < rhs.eigen_value;
}

int ChebSpec::N = -1;
int ChebSpec::L = 0;
std::vector<double> ChebSpec::x = {};
std::vector<std::vector<double> > ChebSpec::Tm = {}, ChebSpec::TmInt = {}, ChebSpec::D = {}, ChebSpec::D2 = {};

ChebSpec::ChebSpec(){}

void ChebSpec::setPotential(const std::vector<double> &XX, const std::vector<double> &VV) 
{
  SolvSpec::setPotential(XX, VV);
  a = XX[0];
  b = XX[XX.size() - 1];
  scal = std::pow((b - a) / 2, 2.0);
}

void ChebSpec::findSpectrum(int nEigen)
{
  if(PotVals.size() < 2) 
  {
    std::cout << "Please use the setPotential() function before using this one." << std::endl;
    return;
  }
  // compute V
  V.resize(L);
  for(int i = 0; i < L; i++)
  {
    V[i] = scal * potFunc.interp(0.5 * ((b - a) * x[i] + b + a));
  }
  // compute Ehat
  Ehat.resize(L);
  for(int i = 0; i < L; i++) 
  {
    Ehat[i].resize(L);
    for(int j = 0; j < L; j++) 
    {
      Ehat[i][j] = D2[i][j] * V[j];
    }
  }
  // compute UE and US
  UE.resize(L);
  US.resize(L);
  for(int i = 0; i < L; i++) 
  {
    UE[i].resize(L);
    US[i].resize(L);
    for(int j = 0; j < L; j++)
    {
      UE[i][j] = 0.5 * (x[i] + 1) * D2[N][j];
      US[i][j] = 0.5 * (x[i] + 1) * Ehat[N][j];
    }
  }
  // compute the A and B matrices, remove the first/last row/column
  int Nr = N - 1;
  arma::mat A(Nr, Nr);
  arma::mat B(Nr, Nr);

  //lapack routines use the column-major order!
  for(int i = 0; i < Nr; i++) 
  {
    for(int j = 0; j < Nr; j++)
    {
      // define A
      A(i, j) = -US[i + 1][j + 1] + Ehat[i + 1][j + 1];
      if(i == j)   // add the identity matrix
      {
        A(i, j) = A[i + j * Nr] - 1;
      }
      // now B
      B(i, j) = D2[i + 1][j + 1] - UE[i + 1][j + 1];
    }
  }

  arma::mat B1A = arma::inv(B) * A;
  arma::cx_vec cxeigval;
  arma::cx_mat cxeigvec;

  // diagonalize
  arma::eig_gen(cxeigval, cxeigvec, B1A);
  arma::mat eigvec = arma::real(cxeigvec);
  arma::vec eigval = arma::real(cxeigval);
  // Sorting the eigensystem by the order of the eigenvalues
  std::vector<EigenPair> Eigensystem(0);
  EigenPair eigpair;
  for(int i = 0; i < Nr; i++)
  {
    eigpair.eigen_value = eigval(i);
    eigpair.eigen_vector = eigvec.col(i);
    Eigensystem.push_back(eigpair);
  }
  std::sort(Eigensystem.begin(), Eigensystem.end(), compEigenPair);
  // now we just need to put everything in our internal format
  // finally safelly add all the modes found to the spectrum (already in a nice way)
  spectrum.clear();
  spectrum.potential = {X, PotVals};
  for (int i = 0; i < nEigen; i++) 
  {
    // build the wavefunction
    std::vector<Point> wf;
    for(int j = 0; j < Nr; j++)
    {
      wf.push_back(Point(0.5 * ((b - a) * x[j] + b + a), Eigensystem[i].eigen_vector(j) ));
    }
    // normalization loop
    double c = 0;
    for(int j = 0; j < Nr; j++)
    {
      c += D[N][j + 1] * wf[j].y * wf[j].y * 0.5 * (b - a);
    }
    // set all the wavefunctions start growing positive from the left
    double s = 1;
    // get the sign of the derivative, this is important since it may be the case
    // that the routine return the same eigenvector with a different sign, in that case
    // the overall coefficient we would like to fit would have change the sign
    for(int j = 4; j < Nr; j++)
    {
      // filter the data and compute the derivative
      double der = (25./12)*int(1e3 * wf[j].y)-4*int(1e3 * wf[j-1].y)+3*int(1e3 * wf[j-2].y)-(4./3) * int(1e3 * wf[j-3].y)+(1./4)* int(1e3 * wf[j-4].y);
      der /= 1e3;
      if(std::fabs(der) > 1e-2) 
      {
        s = der / std::fabs(der);
        break;
      }
    }
    for(int j = 0; j < Nr; j++)
    {
      wf[j].y = s * wf[j].y / std::sqrt(c);
    }
    Mode m(Eigensystem[i].eigen_value / scal, wf);
    spectrum.addMode(m);
  }
}

void ChebSpec::showMatrix(double* A, int Nr)
{
  std::cout << std::endl;
  for(int i = 0; i < Nr; i++) 
  {
    for(int j = 0; j < Nr; j++) 
    {
      std::cout << A[i + j * Nr] << " ";
    }
    std::cout << std::endl;
  }
}

ChebSpec::~ChebSpec(){}

void chebSetN(int n) 
{
  if(ChebSpec::N == n)
    return;
  ChebSpec::N = n;
  ChebSpec::L = n + 1;
  std::cout << "computing chebyshev matrices, N = " << ChebSpec::N << std::endl;
  // initialize x
  ChebSpec::x.clear();
  ChebSpec::x.resize(ChebSpec::L);
  for(int i = 0; i < ChebSpec::L; i++)
  {
    ChebSpec::x[i] = -std::cos(M_PI * i / ChebSpec::N);
  }
  // Tm, the matrix of the Chebyshev polynomial values in the grid. Checked.
  ChebSpec::Tm.clear();
  ChebSpec::Tm.resize(ChebSpec::L);
  for(int i = 0; i < ChebSpec::L; i++) 
  {
    ChebSpec::Tm[i].resize(ChebSpec::L);
    for(int j = 0; j < ChebSpec::L; j++) 
    {
      if(i == 0)      // first
      {  
        ChebSpec::Tm[i][j] = 1;
      }
      else if(i == 1) // second
      {
        ChebSpec::Tm[i][j] = ChebSpec::x[j];
      }
      else
      {
        ChebSpec::Tm[i][j] = 2 * ChebSpec::x[j] * ChebSpec::Tm[i - 1][j] - ChebSpec::Tm[i - 2][j];
      }
    }
  }

  // TmIntInt, the matrix of the integrated from -1 Chebyshev polynomial values in the grid. Checked.
  ChebSpec::TmInt.clear();
  ChebSpec::TmInt.resize(ChebSpec::L);
  for(int i = 0; i < ChebSpec::L; i++) 
  {
    ChebSpec::TmInt[i].resize(ChebSpec::L);
    for(int j = 0; j < ChebSpec::L; j++)
    {
      if(i == 0)      // first
      {
        ChebSpec::TmInt[i][j] = ChebSpec::x[j] + 1;
      }
      else if(i == 1) // second
      {
        ChebSpec::TmInt[i][j] =  0.5 * (ChebSpec::x[j] * ChebSpec::x[j] - 1);
      }
      else if(i < ChebSpec::L - 1)
      {
        ChebSpec::TmInt[i][j] = ChebSpec::Tm[i + 1][j] / (2. * (i + 1)) - ChebSpec::Tm[i - 1][j] / (2. * (i - 1)) + std::pow(-1, i + 1.) / (std::pow(i, 2.) - 1);
      }
      else   // case i = L - 1, right extreme
      {
        ChebSpec::TmInt[i][j] = (2 * ChebSpec::x[j] * ChebSpec::Tm[i][j] - ChebSpec::Tm[i - 1][j]) / (2.0 * (i + 1)) - ChebSpec::Tm[i - 1][j] / (2.0 * (i - 1)) + std::pow(-1, i + 1) / (std::pow(i, 2) - 1);
      }
    }
  }

  // D matrix: the integral from -1 to x. Checked.
  ChebSpec::D.clear();
  ChebSpec::D.resize(ChebSpec::L);
  for(int i = 0; i < ChebSpec::L; i++) 
  {
    ChebSpec::D[i].resize(ChebSpec::L);
    for(int j = 0; j < ChebSpec::L; j++) 
    {
      double s = 0;
      for(int k = 0; k < ChebSpec::L; k++) 
      {
        double f = 1;
        if(k == 0 || k == ChebSpec::L - 1) // prime '' sum
        {
          f = 0.5;
        }
        s += f * ChebSpec::Tm[k][j] * ChebSpec::TmInt[k][i];
      }

      ChebSpec::D[i][j] = (2. / ChebSpec::N) * s;
      // prime '' sum
      if(j == 0 || j == ChebSpec::L - 1)
        ChebSpec::D[i][j] = 0.5 * ChebSpec::D[i][j];
    }
  }
  // finally the D2 = D*D matrix
  // check, the sum of the last row has to be 2 and 0 the sum of the first
  ChebSpec::D2.clear();
  ChebSpec::D2.resize(ChebSpec::L);
  for(int i = 0; i < ChebSpec::L; i++)
  {
    ChebSpec::D2[i].resize(ChebSpec::L);
    for(int j = 0; j < ChebSpec::L; j++)
    {
      double s = 0;
      for(int k = 0; k < ChebSpec::L; k++)
      {
        s += ChebSpec::D[i][k] * ChebSpec::D[k][j];
      }
      ChebSpec::D2[i][j] = s;
    }
  }

  // check, the sum of the last row has to be 2
  double s = 0;
  for(int k = 0; k < ChebSpec::L; k++)
  {
    s += ChebSpec::D2[0][k];
  }
}
