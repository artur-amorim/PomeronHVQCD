#include <iostream>
#include <string>
#include <vector>
#include "HolographicVQCD.h"
#include "HQCDP.h"
#include "GluonKernel.h"
#include "MesonKernel.h"
#include "SigmaProtonProton.h"
#include "methods/optimization/NelderMead.hpp"
#include "schrodinger/chebspec.h"

using namespace std;

int main(int argc, char ** argv)
{
    // Kernel parameters
    double coeff_g, coeff_m;
    int coeff_index;
    // Couplings
    double Im_g0, Im_g1, Im_g2, Im_m0;
    if(argc < 8)
    {
        cout << "Program usage: " << argv[0] << " coeff_g coeff_m coeff_index Im_g0 Im_g1 Im_g2 Im_m0" << endl;
        return 0;
    }
    else
    {
        coeff_g = stod(argv[1]); coeff_m = stod(argv[2]); coeff_index = stoi(argv[3]);
        Im_g0 = stod(argv[4]); Im_g1 = stod(argv[5]); Im_g2 = stod(argv[6]); Im_m0 = stod(argv[7]);
    }
    cout << "Predicting sigma(p p -> X) with" << endl;
    cout << "coeff_index: " << coeff_index << endl;
    cout << "coeff_g: " << coeff_g << " coeff_m: " << coeff_m << endl;
    cout << "Im_g0: " << Im_g0 << " Im_g1: " << Im_g1 << " Im_g2: " << Im_g2 << " Im_m0: " << Im_m0 << endl;

    double mq = hvqcd().QuarkMass();
    SigmaProtonProton sigma("expdata/SigmaProtonProton/SigmaProtonProton_data_W_lt_10000_without_outliers.txt");
    vector<vector<double> > sigma_points = sigma.expKinematics();

    // Setup Gluon and Meson Kernels and GNs vector
    vector<double> gluon_pars = {0, 0, 0, 0, 0, 0, 0};
    vector<double> meson_pars = {0, 0, 0, 0, 0, 0, 0};
    gluon_pars[coeff_index] = coeff_g; meson_pars[coeff_index] = coeff_m;
    GluonKernel gluon(3, gluon_pars);
    MesonKernel meson(1, meson_pars);
    vector<double> GNs = {Im_g0, Im_g1, Im_g2, Im_m0};

    // Setup HQCDP object
    HQCDP hqcdp;
    hqcdp.addProcessObservable(sigma);
    hqcdp.addKernel(gluon);
    hqcdp.addKernel(meson);
    hqcdp.setGNs(GNs);
    // Compute the spectrum
    // initialise Chebyschev matrices
    chebSetN(1000);
    hqcdp.computeSpectrum();
    vector<Spectra> spectrum = hqcdp.getSpectrum();

    // Compute IzNBars
    cout << "Computing sigma(p p -> hadrons) IzNBars" << endl;
    vector<kinStruct> IzNBars = sigma.getIzsBar(sigma_points, spectrum, GNs);
    // Compute IzNs
    cout << "Computing sigma(p p -> hadrons) IzNs" << endl;
    vector<kinStruct> IzNs = sigma.getIzs(sigma_points, spectrum);

    // Compute sigma(p p -> hadrons) chi2
    const int nPoints = sigma_points[0].size();
    double sigma_chi2 = sigma.chi2(IzNs, IzNBars, sigma_points);
    cout << "Number of points: " << nPoints << endl;
    cout << "The sigma(p p -> hadrons) chi2 / points is " << sigma_chi2 / nPoints << endl;
    cout << "The sigma(p p -> hadrons) chi2 / Ndof is " << sigma_chi2 / (nPoints - 5) << endl;

    vector<double> Ws;
    for(double W = 1.5; W < 11000; W += 0.1) Ws.push_back(W);
    vector<double> WPlus(Ws.size(), 0.0), WMinus(Ws.size(), 0.0);
    vector<vector<double>> kinPts = {Ws, WPlus, WMinus};
    // Compute new IzNBars
    IzNBars = sigma.getIzsBar(kinPts, spectrum, GNs);
    // Compute new IzNs
    IzNs = sigma.getIzs(kinPts, spectrum);
    // Compute  predictions of sigma(gamma p -> hadrons)
    std::cout << "Predicting sigma(p p -> hadrons)" << std::endl;
    vector<double> sigma_pred = sigma.predict(IzNs, IzNBars, kinPts, true);

    return 0;
}