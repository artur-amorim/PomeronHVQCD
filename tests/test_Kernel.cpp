#include <iostream>
#include <fstream>
#include <vector>
#include "HolographicVQCD.h"
#include "Reggeon.h"
#include "GluonKernel.h"
#include "schrodinger/schrodinger.h"
#include "schrodinger/chebspec.h"

using namespace std;

int main()
{
    // Compute Chebyschev matrices
    chebSetN(1000);

    // Setup gluon kernel and compute the Reggeons for t = 0
    vector<double> gluon_pars = {1/0.153, 0, 0, 0, 0, 0, 0};
    GluonKernel gluon(4, gluon_pars);

    cout << hvqcd().QuarkMass() << endl;
    /*
    // Ok let's compute the trajectories
    const double jmax = 2.1;
    const double h = 0.025;
    const int n_js = jmax / h ;
    vector<double> js(n_js, 0.0);
    vector<vector<double> > ts(4,std::vector<double>(n_js,0.0));
    for(double i = 0; i < n_js; i++)
    {
        double j = (i+1) * h;
        cout << "Computing potential with j = " << j << endl;
        vector<double> VSch = gluon.computePotential(j);
        vector<double> z = hvqcd().z();
        reverse(begin(z), end(z));
        List  spectrum = computeSpectrum(z, VSch, 4, "cheb");
        js[i] = j;
        for(int k = 0; k < 4; k++) ts[k][i] = spectrum.Es[k];
    }*/

    gluon.computeReggeTrajectories();

    Spline_Interp<double> traj1 = gluon.Trajectory(0);
    Spline_Interp<double> traj2 = gluon.Trajectory(1);
    Spline_Interp<double> traj3 = gluon.Trajectory(2);
    Spline_Interp<double> traj4 = gluon.Trajectory(3);

    cout << "Intercepts:" << endl;
    cout << "j1: " << traj1.interp(0) << endl;
    cout << "j2: " << traj2.interp(0) << endl;
    cout << "j3: " << traj3.interp(0) << endl;
    cout << "j4: " << traj4.interp(0) << endl;
    
    return 0;
}