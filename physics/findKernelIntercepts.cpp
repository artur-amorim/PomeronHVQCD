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
    double coeff_g, coeff_m;
    int g_index, m_index;
    if (argc < 5)
    {
        coeff_g = 10.6221; coeff_m = 5.58333;
        g_index = 2; m_index = 2;
    }
    else
    {
        coeff_g = stod(argv[1]); coeff_m = stod(argv[2]);
        g_index = stoi(argv[3]); m_index = stoi(argv[4]);
    }
    cout << "Starting to hunt for the intercepts with:" << endl;
    cout << "g_index: " << g_index << " m_index: " << m_index << endl;
    cout << "coeff_g: " << coeff_g << " coeff_m: " << coeff_m << endl;

    double mq = hvqcd().QuarkMass();
    vector<double> u = hvqcd().u();
    cout << "Quark mass: " << mq << endl;
    // Compute Chebyschev matrices
    chebSetN(400);

    // Setup gluon kernel
    vector<double> gluon_pars = {0, 0, 0, 0, 0, 0, 0};
    vector<double> meson_pars = {0, 0, 0, 0, 0, 0, 0};
    gluon_pars[g_index] = coeff_g; meson_pars[m_index] = coeff_m;
    GluonKernel gluon(2, gluon_pars);
    MesonKernel meson(1, meson_pars);

    function<double(vector<double>) > func = [&gluon, &meson, &gluon_pars, &meson_pars, g_index, m_index, &u] (const vector<double> &x)
    {
        // Compute Regge Trajectories
        cout << "coeff_g: " << x[0] << " coeff_m: " << x[1] << endl;
        gluon_pars[g_index] = x[0]; meson_pars[m_index] = x[1];
        gluon.computeReggeTrajectories(gluon_pars);
        meson.setKernelPars(meson_pars);
        //meson.computeReggeTrajectories(meson_pars);
        // Compute the intercepts
        Spline_Interp<double> gluon_j0 = gluon.Trajectory(0);
        Spline_Interp<double> gluon_j1 = gluon.Trajectory(1);
        //Spline_Interp<double> meson_j0 = meson.Trajectory(0);
        //Spline_Interp<double> meson_j1 = meson.Trajectory(1);
        double g_j0 = gluon_j0.interp(0), g_j1 = gluon_j1.interp(0); // m_j0 = meson_j0.interp(0);

        // Now we compute the mass of f2
        vector<double> V = meson.computePotential(2);
        const double m_f2 = sqrt(computeSpectrum(u, V, 1).Es[0]);
        // Compute the mass of rho3
        V = meson.computePotential(3);
        const double m_rho3 = sqrt(computeSpectrum(u, V, 1).Es[0]);
        // Compute the mass of f4
        V = meson.computePotential(4);
        const double m_f4 = sqrt(computeSpectrum(u, V, 1).Es[0]);


        // Compute the error squared and return it
        double erms2 = std::fabs(g_j1 - 1.0808) + std::fabs(m_f2 - 1.2755) + std::fabs(m_rho3 - 1.6888) + std::fabs(m_f4 - 2.018);
        cout << "g_j0: " << g_j0 << " g_j1: " << g_j1 << " m_f2: " << m_f2 << " m_rho3: " << m_rho3 << " m_f4: " << m_f4 << endl;
        if (std::isnan(erms2)) erms2 = 1e99;
        cout << "erms2: " << erms2 << endl; 
        return erms2;
    };

    vector<double> x = {coeff_g, coeff_m};
    vector<double> deltas(x.size(), 1);
    x = optimFunction(x, func, deltas, 1e-3);

    std::cout << "Best parameters found for:" << std::endl;
    for(int i = 0; i < x.size(); i++) cout << x[i] << '\t';
    cout << endl;
    
    gluon_pars[g_index] = x[0]; meson_pars[m_index] = x[1];
    gluon.computeReggeTrajectories(gluon_pars);
    meson.computeReggeTrajectories(meson_pars);

    // Compute the intercepts
    Spline_Interp<double> gluon_j0 = gluon.Trajectory(0);
    Spline_Interp<double> gluon_j1 = gluon.Trajectory(1);
    Spline_Interp<double> meson_j0 = meson.Trajectory(0);
    //Spline_Interp<double> meson_j1 = meson.Trajectory(1);
    double g_j0 = gluon_j0.interp(0), g_j1 = gluon_j1.interp(0), m_j0 = meson_j0.interp(0);
    cout << "g_j0: " << g_j0 << " g_j1: " << g_j1 << " m_j0: " << m_j0 << endl;

    // Now we compute the mass of f2
    vector<double> V = meson.computePotential(2);
    const double m_f2 = sqrt(computeSpectrum(u, V, 1).Es[0]);
    // Compute the mass of rho3
    V = meson.computePotential(3);
    const double m_rho3 = sqrt(computeSpectrum(u, V, 1).Es[0]);
    // Compute the mass of f4
    V = meson.computePotential(4);
    const double m_f4 = sqrt(computeSpectrum(u, V, 1).Es[0]);
    cout << "m_f2: " << m_f2 << " m_rho3: " << m_rho3 << " m_f4: " << m_f4 << endl;

    return 0;
}