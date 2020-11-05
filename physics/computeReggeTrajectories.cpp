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

const vector<double> f2Masses = {1.2755, 1.5174, 1.936, 2.011, 2.297, 2.345 };
const vector<double> a2Masses = {1.3169, 1.705};
const vector<double> rho3Masses = {1.6888};
const vector<double> omega3Masses = {1.667};
const vector<double> f4Masses = {2.018};
const vector<double> a4Masses = {1.967};


int main(int argc, char ** argv)
{
    double invls, ag, bg, cg, dg, eg, fg;
    double am, bm, cm, dm, em, fm;
    int N;
    string regge_traj_path;
    if (argc < 16)
    {
        invls = 0; ag = 0; bg = 10.6221; cg = 0; dg = 0; eg = 0; fg = 0;
        am = 0; bm = 5.58333; cm = 0; dm = 0; em = 0; fm = 0;
        N = 400;
        regge_traj_path = "regge_trajectories";
    }
    else
    {
        invls = stod(argv[1]); ag = stod(argv[2]); bg = stod(argv[3]); cg = stod(argv[4]); dg = stod(argv[5]);
        eg = stod(argv[6]); fg = stod(argv[7]);
        am = stod(argv[8]); bm = stod(argv[9]); cm = stod(argv[10]); dm = stod(argv[11]);
        em = stod(argv[12]); fm = stod(argv[13]);
        N = stod(argv[14]); regge_traj_path = argv[15];
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
    chebSetN(N);

    cout << "Computing the Masses of the J = 2 mesons." << endl;
    vector<double> u = hvqcd().u(), Vmeson = meson.computePotential(2);
    std::reverse(u.begin(), u.end());
    List spec = computeSpectrum(u, Vmeson, 6);
    vector<double> m2s = spec.Es;
    cout << "Masses of the tensor mesons f2" << endl;
    for(int i = 0; i < m2s.size(); i++) cout << sqrt(m2s[i]) << '\t';
    cout << endl;

    cout << "Computing the Masses of the J = 3 mesons." << endl;
    Vmeson = meson.computePotential(3);
    spec = computeSpectrum(u, Vmeson, 6);
    m2s = spec.Es;
    cout << "Masses of the rho3, omega3" << endl;
    for(int i = 0; i < m2s.size(); i++) cout << sqrt(m2s[i]) << '\t';
    cout << endl;

    cout << "Computing the Masses of the J = 4 mesons." << endl;
    Vmeson = meson.computePotential(4);
    spec = computeSpectrum(u, Vmeson, 6);
    m2s = spec.Es;
    cout << "Masses of the mesons f4 and a4" << endl;
    for(int i = 0; i < m2s.size(); i++) cout << sqrt(m2s[i]) << '\t';
    cout << endl;

    // Compute the trajectories 
    cout << "Computing the trajectories." << endl;
    Vmeson.clear();
    vector<double> z = hvqcd().z(), Vgluon(0);
    std::reverse(z.begin(), z.end());
    vector<double> ts(0);
    
    ofstream gluon_traj, meson_traj;
    gluon_traj.open(regge_traj_path + "_gluon.txt");
    meson_traj.open(regge_traj_path + "_meson.txt");
    gluon_traj << "j\tt1(j)\tt2(j)\tt3(j)\tt4(j)" << endl;
    meson_traj << "j\tt1(j)\tt2(j)\tt3(j)\tt4(j)" << endl;
    for(double j = 0.1; j <= 4.1; j += 0.01)
    {
        // Compute the potentials
        Vgluon = gluon.computePotential(j);
        Vmeson = meson.computePotential(j);

        // Compute the spectrum of the gluon
        spec = computeSpectrum(z, Vgluon, 4);
        ts = spec.Es;
        gluon_traj << j << '\t' << ts[0] << '\t' << ts[1] << '\t' << ts[2] << '\t' << ts[3] << endl;
        // Compute the spectrum of the meson
        spec = computeSpectrum(u, Vmeson, 4);
        ts = spec.Es;
        meson_traj << j << '\t' << ts[0] << '\t' << ts[1] << '\t' << ts[2] << '\t' << ts[3] << endl;
        
    }
    gluon_traj.close();
    meson_traj.close();

    return 0;
}