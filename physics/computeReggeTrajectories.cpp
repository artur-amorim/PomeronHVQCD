#include <iostream>
#include <fstream>
#include <cmath>
#include "GluonKernel.h"
#include "MesonKernel.h"
#include "HolographicVQCD.h"
#include "schrodinger/chebspec.h"
#include "schrodinger/schrodinger.h"
#include "methods/interpolation/Spline_Interp.hpp"

using namespace std;


int main(int argc, char ** argv)
{
    double invls, ag, bg, cg, dg, eg, fg;
    double am, bm, cm, dm, em, fm;
    int N;
    if (argc < 15)
    {
        invls = 0; ag = 0; bg = 10.6221; cg = 0; dg = 0; eg = 0; fg = 0;
        am = 0; bm = 5.58333; cm = 0; dm = 0; em = 0; fm = 0;
        N = 400;
    }
    else
    {
        invls = stod(argv[1]); ag = stod(argv[2]); bg = stod(argv[3]); cg = stod(argv[4]); dg = stod(argv[5]);
        eg = stod(argv[6]); fg = stod(argv[7]);
        am = stod(argv[8]); bm = stod(argv[9]); cm = stod(argv[10]); dm = stod(argv[11]);
        em = stod(argv[12]); fm = stod(argv[13]);
        N = stod(argv[14]);
    }

    cout << "Computing the gluon Regge trajectories with:" << endl;
    cout << "invls: " << invls << " ag: " << ag << " bg: " << bg << " cg: " << cg << " dg: " << dg;
    cout << " eg: " << eg <<  " fg: " << fg << endl;
    cout << "am: " << am << " bm: " << bm << " cm: " << cm << " dm: " << dm;
    cout << " em: " << em <<  " fm: " << fm << endl;

    vector<double> gluon_pars = {invls, ag, bg, cg, dg, eg, fg};
    GluonKernel gluon(4, gluon_pars);
    vector<double> meson_pars = {invls, am, bm, cm, dm, em, fm};
    MesonKernel meson(2, meson_pars);

    // initialise Chebyschev matrices
    chebSetN(800);

    gluon.computeReggeTrajectories();
    meson.computeReggeTrajectories();

    cout << "Computed the trajectories." << endl;
    // Get the gluon trajectories
    Spline_Interp<double> gluon_traj_0 = gluon.Trajectory(0);
    Spline_Interp<double> gluon_traj_1 = gluon.Trajectory(1);
    Spline_Interp<double> gluon_traj_2 = gluon.Trajectory(2);
    Spline_Interp<double> gluon_traj_3 = gluon.Trajectory(3);
    // Get the meson trajectories
    Spline_Interp<double> meson_traj_0 = meson.Trajectory(0);
    Spline_Interp<double> meson_traj_1 = meson.Trajectory(1);

    cout << "Saving the values in the file regge_trajectories.txt" << endl;
    ofstream myfile;
    myfile.open("regge_trajectories.txt");
    myfile << "t\tgluon_j0\tgluon_j1\tgluon_j2\tgluon_j3\tmeson_j0\tmeson_j1" << endl;
    for(double t = -1; t <= 20; t += 0.01)
    {
        myfile << t << '\t' << gluon_traj_0.interp(t) << '\t' << gluon_traj_1.interp(t) << '\t' << gluon_traj_2.interp(t) << '\t';
        myfile << gluon_traj_3.interp(t) << '\t' << meson_traj_0.interp(t) << '\t' << meson_traj_1.interp(t) << endl; 
    }
    myfile.close();

    cout << "Computation of the Regge Trajectories completed." << endl;

    vector<double> u = hvqcd().u(), Vf2 = meson.computePotential(2);
    std::reverse(u.begin(), u.end());
    List spec = computeSpectrum(u, Vf2, 4);
    vector<double> m2s = spec.Es;
    cout << "Masses of the tensor mesons f2" << endl;
    for(int i = 0; i < m2s.size(); i++) cout << sqrt(m2s[i]) << '\t';
    cout << endl;



    return 0;
}