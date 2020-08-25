#include <iostream>
#include <functional>
#include <cmath>
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
    double ag, am;
    double bg, bm;
    double cg, cm;
    double dg, dm;
    double eg, em;
    double fg, fm;
    if (argc < 3)
    {
        ag = 0; am = 0;
    }
    else
    {
        ag = stod(argv[1]); am = stod(argv[2]);
    }
    cout << "Starting to hunt for the intercepts with:" << endl;
    cout << "ag: " << ag << " am: " << am << endl;

    double mq = hvqcd().QuarkMass();
    // Compute Chebyschev matrices
    chebSetN(400);

    // Setup gluon kernel
    vector<double> gluon_pars = {0, ag, 0, 0, 0, 0, 0};
    vector<double> meson_pars = {0, am, 0, 0, 0, 0, 0};
    GluonKernel gluon(2, gluon_pars);
    MesonKernel meson(1, meson_pars);

    function<double(vector<double>) > func = [&gluon, &meson] (const vector<double> &x){
        // Compute Regge Trajectories
        cout << "ag: " << x[0] << " am: " << x[1] << " eg: " << x[2] << endl;
        gluon.computeReggeTrajectories({0, x[0], 0, 0, 0, 0, 0});
        meson.computeReggeTrajectories({0, x[1], 0, 0, 0, 0, 0});
        // Compute the intercepts
        Spline_Interp<double> gluon_j0 = gluon.Trajectory(0);
        Spline_Interp<double> gluon_j1 = gluon.Trajectory(1);
        Spline_Interp<double> meson_j0 = meson.Trajectory(0);
        //Spline_Interp<double> meson_j1 = meson.Trajectory(1);
        double g_j0 = gluon_j0.interp(0), g_j1 = gluon_j1.interp(0), m_j0 = meson_j0.interp(0);
        // Compute the error squared and return it
        double erms2 = std::fabs(g_j1 - 1.0808) + std::fabs(m_j0 - 0.5475);
        cout << "g_j0: " << g_j0 << " g_j1: " << g_j1 << " m_j0: " << m_j0 << endl;
        cout << "erms2: " << erms2 << endl; 
        return erms2;
    };

    vector<double> x = {ag, am};
    vector<double> deltas(x.size(), 10);
    x = optimFunction(x, func, deltas, 1e-3);

    std::cout << "Best parameters found for:" << std::endl;
    for(int i = 0; i < x.size(); i++) cout << x[i] << '\t';
    cout << endl;
    
    gluon.computeReggeTrajectories({0, x[0], 0, 0, 0, 0, 0});
    meson.computeReggeTrajectories({0, x[1], 0, 0, 0, 0, 0});
    // Compute the intercepts
    Spline_Interp<double> gluon_j0 = gluon.Trajectory(0);
    Spline_Interp<double> gluon_j1 = gluon.Trajectory(1);
    Spline_Interp<double> meson_j0 = meson.Trajectory(0);
    //Spline_Interp<double> meson_j1 = meson.Trajectory(1);
    double g_j0 = gluon_j0.interp(0), g_j1 = gluon_j1.interp(0), m_j0 = meson_j0.interp(0);
    cout << "g_j0: " << g_j0 << " g_j1: " << g_j1 << " m_j0: " << m_j0 << endl;

    return 0;
}