#ifndef HVQCD_H
#define HVQCD_H

#include <boost/numeric/ublas/matrix.hpp>
#include "background.h"

typedef boost::numeric::ublas::matrix<double> matrix_type;

class HVQCD : public Background
{
    private:
        // Parameters of the Vf and k potentials
        double ksc, kU1, W0, WIR, kIR, W1, k1, a1, a2, xf, tau0;
        // Parameters of the w potential
        double wsc, wU1, w0, wIR, w1;
        // Parameters of the Z(lambda) function
        double Za, ca;
        // Quark mass value
        double mq;
        // Flag that says if the background has been computed or not
        bool solved;
        // Containers of the values of tau, dtau/dA, d2q/dA2, d2tau/dA2, d3tau/dA3 and u
        std::vector<double> taus, dtaus, d2qs, d2taus, d3taus, us;
        // Containers with the values of the background fields in the YM regime
        std::vector<double> AYM2, zYM2, uYM2, qYM2, PhiYM2, tauYM2, dqYM2, dPhiYM2, dtauYM2, d2qYM2, d2PhiYM2, d2tauYM2, d3tauYM2;
        // Containers of the values of lambdaUV, taun = exp(A) tau, dlambdaUV/dA, dtaun/dA and AUV
        std::vector<double> lUVs, tauns, dlUVs, dtauns, AUVs;
        // Containers of Astring, dAstring, d2Astring
        std::vector<double> Astrings, dAstrings, d2Astrings;
        // Containers that will be useful latter for Regge theory
        std::vector<double> U2s, aFs, bFs, cFs, dFs, eFs;
        // Containers which contain e^(2 A), e^(2Astring) and l1_2
        std::vector<double> e2As, e2Astrings, l1_2s;
        // Declaration of Vf0 as a function of lambda
        double Vf0l(const double l) const;
        // Declaration of dVf0/dlambda
        double dVf0dlambda(const double l) const;
        // Declaration of d2Vf0/dlambda2
        double d2Vf0dlambda2(const double l) const;
        // Declaration of the Vf potential as a funciton of lambda
        double Vfl(const double l, const double tau) const;
        // Declaration of dVf/dlamda
        double dVfldlambda(const double l, const double tau) const;
        // Declaration of the Vf = Vf0 (1 + tausc tau^2) exp(-tau^2) potential as a funtion of Phi
        double Vf(const double phi, const double tau) const;
        // Declaration of the dVfdPhi = dVf0dPhi exp(-tau^2) potential
        double dVfdPhi(const double phi, const double tau) const;
        // Declaration of the dVfdtau = Vf0 (-2tau)exp(-tau^2) potential
        double dVfdtau(const double phi, const double tau) const;
        // Declaration of the d2VfdPhi2 potential
        double d2VfdPhi2(const double phi, const double tau) const;
        // Declaration of the d2VfdPhidtau = dVf0dPhi (-2tau)exp(-tau^2) potential
        double d2VfdPhidtau(const double phi, const double tau) const;
        // Declaration of the k potential as a function of lambda
        double klambda(const double l) const;
        // Declaration of dk/dlambda potential
        double dkdlambda(const double l) const;
        // Declaration of d2k/dlambda2
        double d2kdlambda2(const double l) const;
        // Declaration of the dw/dlambda potential
        double dwdlambda(const double l) const;
        // Declaration of du/dA
        double dudA(const double q, const double phi, const double dtau, const double A) const;
        // Declaration of dtau/dA when tau > tcut
        double dtauYangMills1(const double q, const double phi, const double tau, const double dphi);
        // Declaration of d2tau/dA2 when tau > tcut
        double d2tauYangMills1(const double q, const double phi, const double tau, const double dq, const double dphi, const double d2phi);
        // Declaration of d2tau/dA2 in the IR
        double d2tauYM2dA2(const double q, const double phi, const double tau, const double dq, const double dphi, const double dtau);
        // Declaration of d3tau/dA3 in the IR
        double d3tauYM(const double q, const double phi, const double tau,
                       const double dq, const double dphi, const double dtau,
                       const double d2q, const double d2phi, const double d2tau);
        // Declaration of the EOMs when tau > tcut
        void eomYangMills1(const state &X, state &dXdA, const double A);
        // Declaration of the EOMs in the IR 
        void eomYangMills2(const state &X, state &dXdA, const double A);
        // Declaration of the observerCoupled function
        void observerYangMills2(const state &X , const double A);
        // Declaration of d2tau/dA2 in the coupled regime
        double d2tauCoupled(const double q, const double phi,const double tau, const double dphi, const double dtau);
        // Declaration of d3tau/dA3 in the coupled regime
        double d3tauCoupled(const double q, const double phi, const double tau, const double dq, const double dphi,
                            const double dtau, const double d2phi, const double d2tau);
        // Declaration of dq/dA in the coupled regime
        double dqCoupled(const double q, const double phi, const double tau, const double dphi, const double dtau);
        // Declaration of d2q/dA2 in the coupled regime
        double d2qCoupled(const double q, const double phi, const double tau,  const double dq, const double dphi,
                          const double dtau, const double d2tau);
        // Declaration of d2Phi/dA2 in the coupled regime
        double d2PhiCoupled(const double q, const double phi, const double tau, const double dphi, const double dtau);
        // Declaration of the jacobian function
        void jacobian(const state &X , matrix_type &jac , const double A, state &dfdt);
        // Declaration of the system of coupled EOMs
        void eomCoupled(const state &X , state &dXdA , const double A);
        // Declaration of the observerCoupled function
        void observerCoupled(const state &X , const double A);
        // Declaration of dq/dA in the UV regime
        double dqUV(const double q, const double lambda , const double dlambda);
        // Declaration of d2q/dA2 in the UV regime
        double d2qUV(const double q, const double lambda , const double dq, const double dlambda, const double d2lambda);
        // Declaration of d2lambda/dA2 in the UV regime
        double d2lambdaUV(const double q, const double lambda , const double dq, const double dlambda);
        // Declaration of d2taun/dA2 in the UV regime
        double d2taunUV(const double q, const double lambda, const double tau, const double dq, const double dlambda, const double dtau, const double A);
        // Declaration of d3taun/dA3 in the UV regime
        double d3taunUV(const double q, const double lambda, const double tau, const double dq, const double dlambda, const double dtau,
                        const double d2q, const double d2lambda, const double d2tau, const double A);
        // Declaration of the EOMs in the UV
        void eomUV(const state &X, state &dXdA, const double A);
        // Declaration of the observerUV function
        void observerUV(const state &X , const double A);
        // Function that processes the final background
        void finalizeBackground();
    public:
        // Class Constructor
        HVQCD(const double ssc = 3.0, const double kksc = 3.0, const double wwsc = 1.56,
             const double WW0 = 2.5, const double ww0 = 1.26, const double kkU1 = 11./9,
             const double wwU1 = 0.0, const double VVgIR = 2.05, const double WWIR = 0.9,
             const double kkIR = 1.8, const double wwIR = 5.0, const double WW1 = 0.0,
             const double kk1 = -0.23, const double ww1 = 0.0, const double aa1 = 0,
             const double aa2 = 1.0, const double xxf = 1.0, const double ttau0 = 1.0, 
             const double za = 133, const double c = 0.26);
        // Class copy constructor
        HVQCD(const HVQCD &hvqcd);
        // Class destructor
        ~HVQCD();
        // Declaration of the Vf0 potential as a function of Phi
        double Vf0(const double phi) const;
        // Declaration of dVf0/dPhi potential
        double dVf0dPhi(const double phi) const;
        // Declaration of the d2Vf0dPhi2 potential
        double d2Vf0dPhi2(const double phi) const;
        // Declaration of the Vtau potential = (1+ tausc tau^2)e^(-tau^2);
        double Vtau(const double tau) const;
        // Declaration of the dVtau potential;
        double dVtau(const double tau) const;
        // Declaration of the Vf = Vf0 Vtau potential as a funtion of Phi
        double Vf(const double phi, const double tau) const;
        // Declaration of the k potential as a function of Phi
        double k(const double phi) const;
        // Declaration of dk/dPhi potential
        double dkdPhi(const double phi) const;
        // Declaration of d2k/dPhi2 potential
        double d2kdPhi2(const double phi) const;
        // Declaration of w(Phi)
        double w(const double phi) const;
        // Declaration of dw/dPhi
        double dwdPhi(const double phi) const;
        // Declaration of the d2w/dPhi2 potential
        double d2wdPhi2(const double phi) const;
        // Declaration of G and its derivatives
        double G(const double q, const double phi, const double dt) const;
        double dG(const double q, const double phi, const double dq, const double dphi, const double dt, const double d2t) const;
        double d2G(const double q, const double phi,
                   const double dq, const double dphi,
                   const double dtau, const double d2q,
                   const double d2phi, const double d2tau,
                   const double d3tau) const;
        // Declaration of Z(lambda)
        double Z(const double l) const;
        // Declaration of the solve function.
        void solve();
        // a1 getter
        double get_a1() const;
        // a2 getter
        double get_a2() const;
        // xf getter
        double get_xf() const;
        // Setter of Za and ca
        void setZa(const double za);
        void setca(const double cca);
        // Declaration of solved getter
        bool isSolved() const;
        // Declaration of Astring getter
        std::vector<double> Astring() const;
        // Declaration of dAstringdz getter
        std::vector<double> dAstring() const;
        // Declaratio of d2Astringdz2 getter
        std::vector<double> d2Astring() const;
        // Declaration of tau getter
        std::vector<double> tau() const;
        // Declaration of dtau/dA getter
        std::vector<double> dtaudA() const;
        // Declaration of d2q/dA2
        std::vector<double> d2qdA2() const;
        // Declaration of d2taudA2
        std::vector<double> d2taudA2() const;
        // Declaration of d3taudA3
        std::vector<double> d3taudA3() const;
        // Declaration of u
        std::vector<double> u() const;
        // Declaration of G
        std::vector<double> G() const;
        // Declaration of U2 getter
        std::vector<double> U2() const;
        // Declaration of aF getter
        std::vector<double> aF() const;
        // Declaration of bF getter
        std::vector<double> bF() const;
        // Declaration of cF getter
        std::vector<double> cF() const;
        // Declaration of dF getter
        std::vector<double> dF() const;
        // Declaration of eF getter
        std::vector<double> eF() const;
        // Declaration of e2A getter
        std::vector<double> e2A() const;
        // Declaration of e2Astring getter
        std::vector<double> e2Astring() const;
        // Declaration of l1_2 = sqrt(exp(Phi))
        std::vector<double> l1_2() const;
        // Declaration of QuarkMass
        double QuarkMass() const;
        // Declaration of saveBackgroundFields
        void saveBackgroundFields(std::string path = "HVQCD_backgrounds.txt");
        // Declaration of savePotentials
        void savePotentials(std::string path = "Potentials.txt");
        // Declaration of the mass squared of the bulk tachyon in the IR
        double TachyonMassSquareIR();
};

// Declaration of the function that computes the potential of Flavour Non-Singlet Vector Mesons
std::vector<double> computeVectorMesonPotential(const HVQCD &hvqcd);
// Declaration of the function that computes the potential of Flavour Non-Singlet Axial Vector Mesons
std::vector<double> computeAxialVectorMesonNonSingletPotential(const HVQCD &hvqcd, const std::vector<double> &VVectorMeson);
// Declaration of the function that computes the potential of Flavour Non-Singlet Pseudoscalar Mesons
std::vector<double> computePseudoScalarMesonPotential(const HVQCD &hvqcd);
// Declaration of the function that computes the potential of Flavour Non-Singlet Scalar Mesons
std::vector<double> computeScalarMesonPotential(const HVQCD &hvqcd);
// Declaration of the function that computes the potential of Flavour Singlet Axial Vector Mesons
std::vector<double> computeAxialVectorMesonSingletPotential(const HVQCD &hvqcd, const std::vector<double> &VAxialVectorMeson);

// Declaration of function that computes the Pseudodoscalar spectrum
std::vector<double> computePseudoScalarMasses(const HVQCD &hvqcd, const int n_masses);

// Save Schrodinger potentials in a given file
void saveSchrodingerPotentials(const HVQCD &hvqcd, std::string path = "SchrodingerPotentials.txt");

// Declaration of computeHVQCDSpectrum function
void computeHVQCDSpectrum(const HVQCD &hvqcd);

// Declaration of computeHVQCDRatios function
void computeHVQCDRatios(const HVQCD &hvqcd);

const std::vector<double> mrhos = {775.5, 1465, 1720, 1909, 2149, 2265};         // Non-Singlet Vector Mesons
const std::vector<double> ma1s = {1230, 1647, 1930, 2096, 2270};                 // Non-Singlet Axial Vector Mesons
const std::vector<double> mpis = {135, 1300, 1812, 2070, 2360};                  // Non-Singlet Pseudoscalar Mesons
const std::vector<double> a0s = {1474, 2025};                                    // Non-Singlet Scalar Mesons
const std::vector<double> mTG = {2150};                                          // Singlet Tensor glueballs
const std::vector<double> momegas = {782.65, 1420, 1670};                        // Singlet Vector Mesons
const std::vector<double> mf1s = {1281.9, 1426.4};                               // Singlet Axial Vector Mesons

// Ratios with the rho vector meson rho
const std::vector<double> RTG_rho = {2.7724};
const std::vector<double> Rrho_rho = {1.8891, 2.2179, 2.4616, 2.7711, 2.9207};
const std::vector<double> Ra1_rho = {1.5861, 2.1238, 2.48872, 2.70277, 2.92714};
const std::vector<double> Rpi_rho = {0.1741,1.6763,2.3366,2.6692, 3.0432};
const std::vector<double> Ra0_rho = {1.9007, 2.6112};
const std::vector<double> Romega_rho = {1.01, 1.83, 2.15};
const std::vector<double> Rf1_rho = {1.65, 1.84};

HVQCD& hvqcd();

#endif