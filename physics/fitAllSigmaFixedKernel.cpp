#include <iostream>
#include <string>
#include <complex>
#include <random>
#include "HolographicVQCD.h"
#include "SigmaGammaProton.h"
#include "SigmaGammaGamma.h"
#include "SigmaProtonProton.h"
#include "GluonKernel.h"
#include "MesonKernel.h"
#include "schrodinger/schrodinger.h"
#include "schrodinger/chebspec.h"
#include "methods/optimization/NelderMead.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    // Gravitational photon couplings
    double k0g, k1g, k0m;
    // Gravitational proton couplings
    double kbar0g, kbar1g, kbar0m;
    string sigma_gp_path, sigma_gg_path, sigma_pp_path;
    if(argc < 10)
    {
        sigma_gg_path = "expdata/SigmaGammaGamma/SigmaGammaGamma_PDG_data_W_gt_4.txt";
        sigma_gp_path = "expdata/SigmaGammaProton/SigmaGammaP_PDG_data_W_gt_461.txt";
        sigma_pp_path = "expdata/SigmaProtonProton/SigmaProtonProton_data_W_lt_10000_without_outliers.txt";
        k0g = 0; k1g = 0; k0m = 0;
        kbar0g = 0; kbar1g = 0; kbar0m = 0;
        cout << "Program usage: " + string(argv[0]) + " sigma_gg_path sigma_gp_path sigma_pp_path k0g k1g k0m kbar0g kbar1g kbar0m " << endl;
        cout << "Using default values." << endl;
    }
    else
    {
        sigma_gg_path = argv[1];
        sigma_gp_path = argv[2]; sigma_pp_path = argv[3];
        k0g = stod(argv[4]); k1g = stod(argv[5]); k0m = stod(argv[6]);
        kbar0g = stod(argv[7]); kbar1g = stod(argv[8]); kbar0m = stod(argv[9]);
    }
    cout << "Starting fit with values:" << endl;
    cout << "k0g: " << k0g << " k1g: " << k1g << " k0m: " << k0m << endl;
    cout << "kbar0g: " << kbar0g << " kbar1g: " << kbar1g << " kbar0m: " << kbar0m << endl;

    // Define the process observables and load the data needed for the fit
    SigmaGammaGamma sigma_gg(sigma_gg_path);
    SigmaGammaProton sigma_gp(sigma_gp_path);
    SigmaProtonProton sigma_pp(sigma_pp_path);

    // Get the experimental points
    int npoints = 0;
    vector<vector<double> > sigma_gg_pts = sigma_gg.expKinematics();
    npoints += sigma_gg_pts[0].size();
    vector<vector<double> > sigma_gp_pts = sigma_gp.expKinematics();
    npoints += sigma_gp_pts[0].size();
    vector<vector<double> > sigma_pp_pts = sigma_pp.expKinematics();
    npoints += sigma_pp_pts[0].size();

    // Setup Chebyschev computation
    chebSetN(400);

    // Setup HardPomeron Kernel and compute the Reggeons
    double bg = -10.6221, bm = -5.58333;
    GluonKernel gluon(2, {0, 0, bg, 0, 0, 0, 0});
    MesonKernel meson(1, {0, 0, bm, 0, 0, 0, 0});
    gluon.computeReggeTrajectories();
    meson.computeReggeTrajectories();
    vector<Reggeon> reggeons_gluon = computeReggeons(gluon, 0, 2);
    vector<Reggeon> reggeons_meson = computeReggeons(meson, 0, 1);
    vector<Reggeon> reggeons = {reggeons_gluon[0], reggeons_gluon[1], reggeons_meson[0]};
    for(int i = 0; i < reggeons.size(); i++) reggeons[i].setIndex(i+1);
    Spectra spec(0, reggeons);
    vector<Spectra> spectrum = {spec};

    // Ok now let's compute all the Izns, because our fitting parameters only enter in the IzNBars
    cout << "Computing sigma(gamma gamma -> X) IzNs" << endl;
    vector<kinStruct> sigma_gg_IzNs = sigma_gg.getIzs(sigma_gg_pts, spectrum);
    cout << "Computing sigma(gamma p -> X) IzNs" << endl;
    vector<kinStruct> sigma_gp_IzNs = sigma_gp.getIzs(sigma_gp_pts, spectrum);
    cout << "Computing sigma(p p -> X) IzNs" << endl;
    vector<kinStruct> sigma_pp_IzNs = sigma_pp.getIzs(sigma_pp_pts, spectrum);

    auto f = [&sigma_gg, &sigma_gp, &sigma_pp, &spectrum, &sigma_gg_pts, &sigma_gp_pts, &sigma_pp_pts,
              &sigma_gg_IzNs, &sigma_gp_IzNs, &sigma_pp_IzNs, &reggeons] (const std::vector<double> & X)
    {
        // X is a vector with the values of k1, k2, k3, k4, kbar1, kbar2, kbar3, kbar4
        vector<double> ks = {X[0], X[1], X[2]};
        vector<double> kbars = {X[3], X[4], X[5]};
        cout << "kg0: " << ks[0] << " kg1: " << ks[1] << " km0: " << ks[2] << " kbarg0: " << kbars[0] << " kbarg1: " << kbars[1] << " kbarm0: " << kbars[2] << endl;
        /*
            Now we need to compute the gn's according to the definition in the notes
        */
        vector<double> Im_gn_gg(reggeons.size()), Im_gn_gp(reggeons.size()), Im_gn_pp(reggeons.size());
        for(int i = 0; i < reggeons.size(); i++)
        {
            double jn = reggeons[i].getJ(), djndt = reggeons[i].getdJdt();
            // Common factor
            complex<double> gn_gg(1/tan(M_PI_2 * jn), 1), gn_gp(1/tan(M_PI_2 * jn), 1), gn_pp(1/tan(M_PI_2 * jn), 1);
            gn_gg = M_PI_2 * gn_gg * djndt / pow(2, jn); 
            gn_gp = M_PI_2 * gn_gp * djndt / pow(2, jn);
            gn_pp = M_PI_2 * gn_pp * djndt / pow(2, jn);
            // Now we make the specific computations
            gn_gg = gn_gg * ks[i]*ks[i] * sigma_gg_IzNs[0].izns[i];
            gn_gp = gn_gp * ks[i] * kbars[i];
            gn_pp = gn_pp * kbars[i] * kbars[i];
            // Take the imaginary part
            Im_gn_gg[i] = imag(gn_gg);
            Im_gn_gp[i] = imag(gn_gp);
            Im_gn_pp[i] = imag(gn_pp);
        }
                   
        // Compute the IzNsBar
        vector<kinStruct> sigma_gg_IzNBars = sigma_gg.getIzsBar(sigma_gg_pts, spectrum, Im_gn_gg);
        vector<kinStruct> sigma_gp_IzNBars = sigma_gp.getIzsBar(sigma_gp_pts, spectrum, Im_gn_gp);
        vector<kinStruct> sigma_pp_IzNBars = sigma_pp.getIzsBar(sigma_pp_pts, spectrum, Im_gn_pp);

        // Compute the chi2 of each process
        double sigma_gg_chi2 = sigma_gg.chi2(sigma_gg_IzNs, sigma_gg_IzNBars, sigma_gg_pts);
        double sigma_gp_chi2 = sigma_gp.chi2(sigma_gp_IzNs, sigma_gp_IzNBars, sigma_gp_pts);
        double sigma_pp_chi2 = sigma_pp.chi2(sigma_pp_IzNs, sigma_pp_IzNBars, sigma_pp_pts);

        cout << "Im_gn_gg:" << endl;
        for(int i = 0; i < Im_gn_gg.size(); i++) cout << Im_gn_gg[i] << '\t';
        cout << endl;
        cout << "Im_gn_gp:" << endl;
        for(int i = 0; i < Im_gn_gp.size(); i++) cout << Im_gn_gp[i] << '\t';
        cout << endl;
        cout << "Im_gn_pp:" << endl;
        for(int i = 0; i < Im_gn_pp.size(); i++) cout << Im_gn_pp[i] << '\t';
        cout << endl;

        double chi2 = sigma_gg_chi2 + sigma_gp_chi2 + sigma_pp_chi2;
        cout << "chi2: " << chi2 << endl; 
        return chi2;
    };


    // Start the fit now
    vector<double> X_guess = {k0g, k1g, k0m, kbar0g, kbar1g, kbar0m};
    vector<double> delta(X_guess.size(), 10);
    vector<double> X_opt = optimFunction(X_guess, f, delta, 1e-12);
    
    // Print X_opt
    cout << "Found a minimum for: ";
    for(int i = 0; i < X_opt.size(); i++) cout << X_opt[i] << " ";
    cout << endl;
    cout << "Number of experimental points: " << npoints << endl;
    cout << "Number of degrees of freedom: " << npoints - X_opt.size() << endl;
    cout << "Best chi2 / Ndof: " << f(X_opt) / (npoints - X_opt.size()) << endl;

    return 0;
}