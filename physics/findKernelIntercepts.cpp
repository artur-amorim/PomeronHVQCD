#include <iostream>
#include <functional>
#include "HolographicVQCD.h"
#include "GluonKernel.h"
#include "MesonKernel.h"
#include "schrodinger/chebspec.h"
#include "schrodinger/schrodinger.h"
#include "methods/optimization/NelderMead.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    double invls;
    if (argc < 2)
    {
        invls = 0.254119;
    }
    else
    {
        invls = stod(argv[1]);
    }
    cout << "Starting to hunt for the intercepts with:" << endl;
    cout << "invls: " << invls << endl;

    double mq = hvqcd().QuarkMass();
    // Compute Chebyschev matrices
    chebSetN(1000);

    // Setup gluon kernel
    vector<double> gluon_pars = {invls, 0, 0, 0, 0, 0, 0};
    GluonKernel gluon(2, gluon_pars);
    MesonKernel meson(1, gluon_pars);

    function<double(vector<double>) > func = [&gluon, &meson] (const vector<double> &x){
        // Compute Regge Trajectories
        gluon.computeReggeTrajectories({x[0], 0, 0, 0, 0, 0, 0});
        meson.computeReggeTrajectories({x[0], 0, 0, 0, 0, 0, 0});
        // Compute the intercepts
        Spline_Interp<double> gluon_j0 = gluon.Trajectory(0);
        Spline_Interp<double> gluon_j1 = gluon.Trajectory(1);
        Spline_Interp<double> meson_j0 = meson.Trajectory(0);
        double g_j0 = gluon_j0.interp(0), g_j1 = gluon_j1.interp(0), m_j0 = meson_j0.interp(0);
        // Compute the error squared and return it
        double erms2 = pow(g_j0 - 1.4, 2) + pow(g_j1 - 1.08, 2) + pow(m_j0 - 0.45, 2);
        cout << "invls: " << x[0] << endl;
        cout << "g_j0: " << g_j0 << " g_j1: " << g_j1 << " m_j0: " << m_j0 << endl;
        cout << "erms2: " << erms2 << endl; 
        return erms2;
    };

    gluon_pars = optimFunction(gluon_pars, func, 10, 1e-3);

    std::cout << "Best parameters found for:" << std::endl;
    std::cout << "invls:\t" << gluon_pars[0] << std::endl;
    
    gluon.computeReggeTrajectories({gluon_pars[0],0,0,0,0,0,0});
    meson.computeReggeTrajectories({gluon_pars[0],0,0,0,0,0,0});
    // Compute the intercepts
    Spline_Interp<double> gluon_j0 = gluon.Trajectory(0);
    Spline_Interp<double> gluon_j1 = gluon.Trajectory(1);
    Spline_Interp<double> meson_j0 = meson.Trajectory(0);

    double g_j0 = gluon_j0.interp(0), g_j1 = gluon_j1.interp(0), m_j0 = meson_j0.interp(0);
    cout << "g_j0: " << g_j0 << "g_j1: " << g_j1 << " m_j0: " << m_j0 << endl;

    return 0;
}