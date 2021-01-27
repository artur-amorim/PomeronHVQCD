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
    double invls_g, ag, bg, cg, dg, eg, fg;
    double invls_m, am, bm, cm, dm, em, fm;
    int N;
    string regge_traj_path;
    if (argc < 17)
    {
        invls_g = 0; ag = 0; bg = 10.6221; cg = 0; dg = 0; eg = 0; fg = 0;
        invls_m = 0; am = 0; bm = 5.58333; cm = 0; dm = 0; em = 0; fm = 0;
        N = 400;
        regge_traj_path = "regge_trajectories";
    }
    else
    {
        invls_g = stod(argv[1]); ag = stod(argv[2]); bg = stod(argv[3]); cg = stod(argv[4]); dg = stod(argv[5]);
        eg = stod(argv[6]); fg = stod(argv[7]);
        invls_m = stod(argv[8]); am = stod(argv[9]); bm = stod(argv[10]); cm = stod(argv[11]); dm = stod(argv[12]);
        em = stod(argv[13]); fm = stod(argv[14]);
        N = stod(argv[15]); regge_traj_path = argv[16];
    }

    cout << "Computing the gluon Regge trajectories with:" << endl;
    cout << "invls_g: " << invls_g << " ag: " << ag << " bg: " << bg << " cg: " << cg << " dg: " << dg;
    cout << " eg: " << eg <<  " fg: " << fg << endl;
    cout << "invls_m: " << invls_m << " am: " << am << " bm: " << bm << " cm: " << cm << " dm: " << dm;
    cout << " em: " << em <<  " fm: " << fm << endl;

    double mq = hvqcd().QuarkMass();
    
    vector<double> gluon_pars = {invls_g, ag, bg, cg, dg, eg, fg};
    GluonKernel gluon(4, gluon_pars);
    vector<double> meson_pars = {invls_m, am, bm, cm, dm, em, fm};
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

    vector<double> js(0);
    vector<vector<double>> ts_gluon(4, vector<double>(0));
    vector<vector<double>> ts_meson(4, vector<double>(0));

    for(double j = 0.1; j <= 4.1; j += 0.01)
    {
        // Append value of j
        js.push_back(j);
        // Compute the potentials
        Vgluon = gluon.computePotential(j);
        Vmeson = meson.computePotential(j);

        // Compute the spectrum of the gluon
        spec = computeSpectrum(z, Vgluon, 4);
        ts = spec.Es;
        gluon_traj << j << '\t' << ts[0] << '\t' << ts[1] << '\t' << ts[2] << '\t' << ts[3] << endl;
        // Append the t values to ts_gluon
        ts_gluon[0].push_back(ts[0]); ts_gluon[1].push_back(ts[1]); 
        ts_gluon[2].push_back(ts[2]); ts_gluon[3].push_back(ts[3]);
        
        // Compute the spectrum of the meson
        spec = computeSpectrum(u, Vmeson, 4);
        ts = spec.Es;
        meson_traj << j << '\t' << ts[0] << '\t' << ts[1] << '\t' << ts[2] << '\t' << ts[3] << endl;
        // Append the t values to ts_meson
        ts_meson[0].push_back(ts[0]); ts_meson[1].push_back(ts[1]); 
        ts_meson[2].push_back(ts[2]); ts_meson[3].push_back(ts[3]);
        
    }
    gluon_traj.close();
    meson_traj.close();

    // Compute the intercepts of the Pomeron and meson trajectories
    Spline_Interp<double> gluon_traj_0(ts_gluon[0], js), gluon_traj_1(ts_gluon[1], js), gluon_traj_2(ts_gluon[2], js), gluon_traj_3(ts_gluon[3], js) ;
    Spline_Interp<double> meson_traj_0(ts_meson[0], js), meson_traj_1(ts_meson[1], js), meson_traj_2(ts_meson[2], js), meson_traj_3(ts_meson[3], js) ;

    const double j0g = gluon_traj_0.interp(0), j1g = gluon_traj_1.interp(0), j2g = gluon_traj_2.interp(0), j3g = gluon_traj_3.interp(0);
    const double j0m = meson_traj_0.interp(0), j1m = meson_traj_1.interp(0), j2m = meson_traj_2.interp(0), j3m = meson_traj_3.interp(0);

    cout << "Gluon Intercepts:" << endl;
    cout << "j0g = " << j0g << " j1g = " << j1g << " j2g = " << j2g << " j3g = " << j3g << endl;
    cout << "Meson Intercepts:" << endl;
    cout << "j0m = " << j0m << " j1m = " << j1m << " j2m = " << j2m << " j3m = " << j3m << endl;

    return 0;
}