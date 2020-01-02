#include <iostream>
#include <vector>
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

    // Print Eigenvalues obtainded by Chebyschev method
    cout << "Eigenvalues obtained using Chebyschev method:" << endl;
    for(int i = 0; i < nEigenvalues; i++) cout << chebEigenVals[i] << '\t';
    cout << endl;

    // Print Eigenvalues obtainded by Numerov method
    cout << "Eigenvalues obtained using Numerov method:" << endl;
    for(int i = 0; i < nEigenvalues; i++) cout << numerovEigenVals[i] << '\t';
    cout << endl;

    return 0;
}