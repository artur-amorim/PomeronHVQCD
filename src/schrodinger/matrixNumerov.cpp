#include <armadillo>
#include "schrodinger/matrixNumerov.h"

MatrixNumerov::MatrixNumerov(): SolvSpec() {}

void MatrixNumerov::findSpectrum(const int nEigen)
{
    // Computing d
    const double d = 0.01 ;

    // Get xmin and xmax from X. We assume it is
    // in ascending order
    const double xmin = X[0];
    const double xmax = X[X.size()-1];
    // Compute N
    const int N = (xmax - xmin) / d - 1;
    // compute V1, V2, ... VN
    arma::mat V_matrix = arma::zeros(N,N);
    for(int i = 0; i < N; i++) V_matrix(i,i) = potFunc.interp(xmin + (i+1) * d);

    // Define A and B matrices
    arma::mat A = arma::zeros(N,N), B = arma::zeros(N,N);
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            if(j == i)
            {
                // Diagonal elements
                A(i, j) = - 2.0 / (d*d);
                B(i, j) = 10.0/12.0;
            }
            else if(j == i - 1 || j == i + 1)
            {
                // Non-zero off-diagonal elements
                A(i,j) = 1.0 / (d*d);
                B(i,j) = 1.0 / 12.0 ;
            }
        }
    }
    // Compute KinMatrix
    const arma::mat KinMatrix = arma::pinv(B) * A;

    // Compute the Hamiltonian matrix H
    const arma::mat H = -KinMatrix + V_matrix;

    // Compute the eigenvalues of H
    arma::vec cxeigval;
    arma::mat cxeigvec;

    arma::eig_sym(cxeigval, cxeigvec, H);
    arma::mat eigvec = arma::real(cxeigvec);
    arma::vec eigval = arma::real(cxeigval);

    // Sorting the eigensystem using ascending order of the eigenvalues
    std::vector<EigenPair> Eigensystem(0);
    EigenPair eigpair;
    for(int i = 0; i < N; i++)
    {
        eigpair.eigen_value = eigval(i);
        eigpair.eigen_vector = eigvec.col(i);
        Eigensystem.push_back(eigpair);
    }
    sort(Eigensystem.begin(), Eigensystem.end(), compEigenPair);

    spectrum.clear();
    spectrum.potential = {X, PotVals};
    for (int i = 0; i < nEigen; i++) 
    {
        // build the wavefunction
        std::vector<Point> wf;
        wf.push_back(Point(xmin, 0));
        for(int j = 0; j < N; j++)
        {
            wf.push_back(Point(xmin + (i + 1) * d, Eigensystem[i].eigen_vector(j) ));
        }
        wf.push_back(Point(xmax, 0));

        // Add mode to the spectrum
        Mode m(Eigensystem[i].eigen_value, wf);
        spectrum.addMode(m);
    }
}

MatrixNumerov::~MatrixNumerov(){}
