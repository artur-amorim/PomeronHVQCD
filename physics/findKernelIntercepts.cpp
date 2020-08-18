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
    GluonKernel gluon(1, gluon_pars);
    MesonKernel meson(1, gluon_pars);

    function<double(vector<double>) > func = [&gluon, &meson] (const vector<double> &params){
        // Compute Regge Trajectories
        gluon.computeReggeTrajectories(params);
        meson.computeReggeTrajectories(params);
        // Compute the intercepts
        Spline_Interp<double> gluon_j0 = gluon.Trajectory(0);
        Spline_Interp<double> meson_j0 = meson.Trajectory(0);
        double g_j0 = gluon_j0.interp(0), m_j0 = meson_j0.interp(0);
        // Compute the error squared and return it
        double erms2 = pow(g_j0 - 1.08, 2) + pow(m_j0 - 0.45, 2);
        cout << "g_j0: " << g_j0 << " m_j0: " << m_j0 << endl;
        cout << "erms2: " << erms2 << endl; 
        return erms2;
    };

    gluon_pars = optimFunction(gluon_pars, func, 10, 1e-3);

    std::cout << "Best parameters found for:" << std::endl;
    std::cout << "invls:\t" << gluon_pars[0] << std::endl;
    std::cout << "a:\t"     << gluon_pars[1] << std::endl;
    std::cout << "b:\t"     << gluon_pars[2] << std::endl;
    std::cout << "c:\t"     << gluon_pars[3] << std::endl;
    std::cout << "d:\t"     << gluon_pars[4] << std::endl;
    std::cout << "e:\t"     << gluon_pars[5] << std::endl;
    std::cout << "f:\t"     << gluon_pars[6] << std::endl;

    gluon.computeReggeTrajectories(gluon_pars);
    meson.computeReggeTrajectories(gluon_pars);
    // Compute the intercepts
    Spline_Interp<double> gluon_j0 = gluon.Trajectory(0);
    Spline_Interp<double> meson_j0 = meson.Trajectory(0);

    double g_j0 = gluon_j0.interp(0), m_j0 = meson_j0.interp(0);

    cout << "g_j0: " << g_j0 << " m_j0: " << m_j0 << endl;

    return 0;
}