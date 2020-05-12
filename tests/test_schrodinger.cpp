#include <iostream>
#include <vector>
#include "schrodinger/chebspec.h"
#include "schrodinger/schrodinger.h"

using namespace std;

int main()
{
    int nEigenvalues = 10;
    
    // Compute the potential
    vector<double> Z, V;
    for(double z = -20; z <= 20; z += 0.1)
    {
        Z.push_back(z);
        V.push_back(z * z);
    }

    // Compute Chebyschev matrices
    chebSetN(800);

    // Compute spectrum using Chebyshev method
    cout << "Computing spectrum using Chebyschev method" << endl;
    List chebSpec = computeSpectrum(Z , V, nEigenvalues, "cheb");
    vector<double> chebEigenVals = chebSpec.Es;

    // Compute spectrum using Numerov method
    cout << "Computing spectrum using Numerov method" << endl;
    List numerovSpec = computeSpectrum(Z, V, nEigenvalues, "numerov");
    vector<double> numerovEigenVals = numerovSpec.Es;

    // Compute spectrum using Matrix Numerov method
    cout << "Computing spectrum using Matrix Numerov method" << endl;
    List matrix_numerov_Spec = computeSpectrum(Z, V, nEigenvalues, "matrix_numerov");
    vector<double> matrix_numerov_EigenVals = matrix_numerov_Spec.Es;

    // Print Eigenvalues obtainded by Chebyschev method
    cout << "Eigenvalues obtained using Chebyschev method:" << endl;
    for(int i = 0; i < nEigenvalues; i++) cout << chebEigenVals[i] << '\t';
    cout << endl;

    // Print Eigenvalues obtainded by Numerov method
    cout << "Eigenvalues obtained using Numerov method:" << endl;
    for(int i = 0; i < nEigenvalues; i++) cout << numerovEigenVals[i] << '\t';
    cout << endl;

    // Print Eigenvalues obtainded by Matrix Numerov method
    cout << "Eigenvalues obtained using Matrix Numerov method:" << endl;
    for(int i = 0; i < nEigenvalues; i++) cout << matrix_numerov_EigenVals[i] << '\t';
    cout << endl;

    return 0;
}