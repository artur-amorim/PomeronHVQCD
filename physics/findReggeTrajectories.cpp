#include <iostream>
#include <functional>
#include "HolographicVQCD.h"
#include "GluonKernel.h"
#include "schrodinger/chebspec.h"
#include "schrodinger/schrodinger.h"
#include "methods/optimization/NelderMead.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    double invls, a, b, c, d, e, f;
    if (argc < 8)
    {
        invls = 0.254119; a = 13.9538; b = 0.921665; c = 2.03904; d = -2.7305; e = -0.473787;
        f = -0.517072;
    }
    else
    {
        invls = stod(argv[1]); a = stod(argv[2]); b = stod(argv[3]); c = stod(argv[4]);
        d = stod(argv[5]); e = stod(argv[6]); f = stod(argv[7]);
    }
    cout << "Starting to hunt the intercepts with:" << endl;
    cout << "invls: " << invls << " a: " << a << " b: " << b << " c: " << c << " d: " << d;
    cout << " e: " << e << " f: " << f << endl;
    double mq = hvqcd().QuarkMass();
    // Compute Chebyschev matrices
    chebSetN(400);

    // Setup gluon kernel
    vector<double> gluon_pars = {invls, a, b, c, d, e, f};
    GluonKernel gluon(4, gluon_pars);

    function<double(vector<double>) > func = [&gluon] (const vector<double> &params){
        // Compute Regge Trajectories
        gluon.computeReggeTrajectories(params);
        // Compute the intercepts
        Spline_Interp<double> traj0 = gluon.Trajectory(0);
        Spline_Interp<double> traj1 = gluon.Trajectory(1);
        Spline_Interp<double> traj2 = gluon.Trajectory(2);
        Spline_Interp<double> traj3 = gluon.Trajectory(3);
        double j0 = traj0.interp(0), j1 = traj1.interp(0);
        double j2 = traj2.interp(0), j3 = traj3.interp(0);
        // Compute the error squared and return it
        double erms2 = pow(j0 - 1.17, 2) + pow(j1 - 1.09, 2) + pow(j2 - 0.969, 2) + pow(j3 - 0.900, 2);
        cout << "j0: " << j0 << " j1: " << j1 << " j2: " << j2 << " j3: " << j3 << endl;
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
    // Compute the intercepts
    Spline_Interp<double> traj0 = gluon.Trajectory(0);
    Spline_Interp<double> traj1 = gluon.Trajectory(1);
    Spline_Interp<double> traj2 = gluon.Trajectory(2);
    Spline_Interp<double> traj3 = gluon.Trajectory(3);
    double j0 = traj0.interp(0), j1 = traj1.interp(0);
    double j2 = traj2.interp(0), j3 = traj3.interp(0);

    cout << "j0: " << j0 << " j1: " << j1 << " j2: " << j2 << " j3: " << j3 << endl;

    return 0;
}