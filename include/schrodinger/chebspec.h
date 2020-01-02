#ifndef CHEBSPEC_H
#define CHEBSPEC_H

#include <armadillo>
#include "solvspec.h"

struct EigenPair
{
    double eigen_value;
    arma::vec eigen_vector;
    EigenPair(const double eigval = 0, const arma::vec & eigvec = arma::vec());
};

bool compEigenPair(const EigenPair &lhs, const EigenPair &rhs);

class ChebSpec : public SolvSpec{
    private:
        double a, b, scal;
        std::vector<double> V;
        std::vector< std::vector<double> > Ehat, UE, US;
        void invMat(double*, int);
        void showMatrix(double* A, int Nr);
    public:
        static int N, L;
        static std::vector<double> x;
        static std::vector< std::vector<double> > Tm, TmInt, D, D2;
        ChebSpec();
        void findSpectrum(int nEigen);
        void setPotential(const std::vector<double> &XX, const std::vector<double> &VV);
        ~ChebSpec();
};

void chebSetN(int n);
#endif