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
    double invls, ag, bg, cg, dg, eg, fg;
    double am, bm, cm, dm, em, fm;
    if (argc < 14)
    {
        invls = -0.0524916; ag = 0.263489; bg = -10.6215; cg = 5.86537; dg = 4.41891; eg = -1.3621; fg = -0.158571;
        am = 1.58592; bm = -5.57058; cm = 10.2666; dm = -0.110224; em = 1.9825; fm = -3.20947;
    }
    else
    {
        invls = stod(argv[1]); ag = stod(argv[2]); bg = stod(argv[3]); cg = stod(argv[4]); dg = stod(argv[5]);
        eg = stod(argv[6]); fg = stod(argv[7]);
        am = stod(argv[8]); bm = stod(argv[9]); cm = stod(argv[10]); dm = stod(argv[11]);
        em = stod(argv[12]); fm = stod(argv[13]);
    }
    cout << "Starting to hunt the intercepts with:" << endl;
    cout << "invls: " << invls << " ag: " << ag << " bg: " << bg << " cg: " << cg << " dg: " << dg;
    cout << " eg: " << eg <<  " fg: " << fg << endl;
    cout << "am: " << am << " bm: " << bm << " cm: " << cm << " dm: " << dm;
    cout << " em: " << em <<  " fm: " << fm << endl;
    double mq = hvqcd().QuarkMass();
    // Compute Chebyschev matrices
    chebSetN(400);

    // Setup gluon and meson kernels
    GluonKernel gluon(2, {invls, ag, bg, cg, dg, eg, fg});
    MesonKernel meson(1, {invls, am, bm, cm, dm, em, fm});

    function<double(vector<double>) > func = [&gluon, &meson] (const vector<double> &x){
        
        // Compute Regge Trajectories
        gluon.computeReggeTrajectories({x[0], x[1], x[2], x[3], x[4], x[5], x[6]});
        meson.computeReggeTrajectories({x[0], x[7], x[8], x[9], x[10], x[11], x[12]});
        // Compute the intercepts
        Spline_Interp<double> traj0 = gluon.Trajectory(0);
        Spline_Interp<double> traj1 = gluon.Trajectory(1);
        Spline_Interp<double> traj2 = meson.Trajectory(0);
        double j1g = traj0.interp(0), j2g = traj1.interp(0);
        double j1m = traj2.interp(0);
        // Compute the error squared and return it
        double erms2 = fabs(j1g - 1.4) + fabs(j2g - 1.0808) + fabs(j1m - 0.5475);
        cout << "j1g: " << j1g << " j2g: " << j2g << " j1m: " << j1m << endl;
        cout << "erms2: " << erms2 << endl; 
        return erms2;
    };

    vector<double> kernels_pars = optimFunction({invls, ag, bg, cg, dg, eg, fg, am, bm, cm, dm, em, fm}, func, 10, 1e-12);

    std::cout << "Best parameters found for:" << std::endl;
    std::cout << "invls:\t" << kernels_pars[0] << std::endl;
    std::cout << "ag:\t"     << kernels_pars[1] << std::endl;
    std::cout << "bg:\t"     << kernels_pars[2] << std::endl;
    std::cout << "cg:\t"     << kernels_pars[3] << std::endl;
    std::cout << "dg:\t"     << kernels_pars[4] << std::endl;
    std::cout << "eg:\t"     << kernels_pars[5] << std::endl;
    std::cout << "fg:\t"     << kernels_pars[6] << std::endl;
    std::cout << "am:\t"     << kernels_pars[7] << std::endl;
    std::cout << "bm:\t"     << kernels_pars[8] << std::endl;
    std::cout << "cm:\t"     << kernels_pars[9] << std::endl;
    std::cout << "dm:\t"     << kernels_pars[10] << std::endl;
    std::cout << "em:\t"     << kernels_pars[11] << std::endl;
    std::cout << "fm:\t"     << kernels_pars[12] << std::endl;

    gluon.computeReggeTrajectories({kernels_pars[0], kernels_pars[1], kernels_pars[2], kernels_pars[3], kernels_pars[4], kernels_pars[5], kernels_pars[6]});
    gluon.computeReggeTrajectories({kernels_pars[0], kernels_pars[7], kernels_pars[8], kernels_pars[9], kernels_pars[10], kernels_pars[11], kernels_pars[12]});

    // Compute the intercepts
    Spline_Interp<double> traj0 = gluon.Trajectory(0);
    Spline_Interp<double> traj1 = gluon.Trajectory(1);
    Spline_Interp<double> traj2 = meson.Trajectory(0);
    double j1g = traj0.interp(0), j2g = traj1.interp(0);
    double j1m = traj2.interp(0);

    cout << "j1g: " << j1g << " j2g: " << j2g << " j1m: " << j1m << endl;

    return 0;
}