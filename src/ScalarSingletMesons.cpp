#ifndef _SCALARSINGLET_CPP
#define _SCALARSINGLET_CPP

#include <cmath>
#include "HolographicVQCD.h"

// Function definitions to compute Scalar Singlet states
double kSS(const double q, const double dq, const double dphi, const double d2phi)
{
    // k(A) function that appears in the coupled EOMs of scalar singlet fluctuations
    double ans = 4.0 - dq/q + 2.0 * d2phi/dphi;
    return ans;
}

double pSS(const double W0, const double xf, const double sc, const double ksc, const double kU1, const double kIR, const double WIR, const double k1, const double W1, const double q, const double phi, const double tau, const double dphi, const double dtau)
{
    // p(A) function that appears in the coupled EOMs of scalar singlet fluctuations
    double vf = Vf(W0, xf, sc, WIR, W1, phi, tau);
    double dvf = dVfdPhi(W0, xf, sc, WIR, W1, phi, tau);
    double kPhi = k(W0, xf, sc, ksc, kU1, kIR, k1, phi);
    double dkPhi = dkdPhi(W0, xf, sc, ksc, kU1, kIR, k1, phi);
    double ans = -vf*(kPhi*kPhi*pow(dtau,4)/(3.0*q*q*(1+kPhi*dtau*dtau/(q*q))) + 3*dkPhi*dtau*dtau*(2+kPhi*dtau*dtau/(q*q))/(8.0*(1+kPhi*dtau*dtau/(q*q))*dphi) + 3.0*kPhi*dtau*dtau*dvf/(4.0*dphi))/(2.0*sqrt(1+kPhi*dtau*dtau/(q*q))) ;
    return ans;
}

double qSS(const double W0, const double xf, const double sc, const double ksc, const double kU1, const double kIR, const double WIR, const double k1, const double W1, const double q, const double phi, const double tau, const double dphi, const double dtau)
{
    // q(A) function that appears in the coupled EOMs of scalar singlet fluctuations
    double vf0 = Vf0(W0, xf, sc, WIR, W1, phi);
    double dvf0 = dVf0dPhi(W0, xf, sc, WIR, W1, phi);
    double kPhi = k(W0, xf, sc, ksc, kU1, kIR, k1, phi);
    double dkPhi = dkdPhi(W0, xf, sc, ksc, kU1, kIR, k1, phi);
    double ans = dvf0*(1+kPhi*dtau*dtau/(q*q))/vf0 + 0.5*dkPhi*(2+kPhi*dtau*dtau/(q*q))/kPhi;
    ans = ans + (4./9)*kPhi*dtau*dtau*dphi/(q*q);
    ans = dphi * ans;
    return ans;
}

double tSS(const double W0, const double xf, const double sc, const double ksc, const double kU1, const double kIR, const double k1, const double q, const double phi, const double dt, const double A)
{
    // t(A) function that appears in the coupled EOMs os scalar singlet perturbations
    double kPhi = k(W0, xf, sc, ksc, kU1, kIR, k1, phi);
    return exp(-2*A) * (q*q + kPhi * dt * dt);
}

double nSS(const double W0, const double xf, const double sc, const double ksc, const double kU1, const double kIR, const double WIR, const double k1, const double W1, const double q, const double phi, const double tau, const double dq, const double dphi, const double dtau)
{
    double kPhi = k(W0, xf, sc, ksc, kU1, kIR, k1, phi);
    double dkPhi = dkdPhi(W0, xf, sc, ksc, kU1, kIR, k1, phi);
    double vf0 = Vf0(W0, xf, sc, WIR, W1, phi);
    double dvf0 = dVf0dPhi(W0, xf, sc, WIR, W1, phi);
    double ans = -4 + dq/q - 4 * q * q * tau / (kPhi * dtau) + 4 * kPhi * pow(dtau/q,2);
    ans = ans - dkPhi * dphi / kPhi - dphi * dvf0 / vf0 + 0.5 * dkPhi * dphi * pow(dtau/q, 2);
    ans = ans + kPhi * dphi * dtau * dtau * dvf0 / (q*q*vf0);
    return ans;
}

double RicciA(const double q, const double dq)
{
    // Given q and dq as funtions of A
    // Returns the value of the Ricci Scalar
    double ans = -20/pow(q,2) + 8 * dq / pow(q,3);
    return ans;
}

double RA(const double W0, const double xf, const double sc, const double ksc, 
          const double kU1, const double kIR, const double k1, const double q, 
          const double phi, const double dtau)
{
    // Returns the R function that appears in the Scalar Singlet EOMs
    // Given q, phi and dtau as functions of A
    double kPhi = k(W0, xf, sc, ksc, kU1, kIR, k1, phi);
    double ans = sqrt(1+kPhi*pow(dtau/q, 2));
    return ans;
}

double QA(const double sc, const double VgIR, const double q,
          const double phi, const double dphi, const double A)
{
    // Returns the Q function that appears in the Scalar Singlet EOMs
    // Given q, phi and dphi as functions of A
    double vg = Vg(sc, VgIR, phi);
    double ans = exp(2*A) * (vg - pow(dphi/q,2));
    return ans;
}

double D1A(const double W0, const double xf, const double sc,
           const double WIR, const double W1, const double ksc,
           const double kU1, const double kIR, const double k1,
           const double q, const double phi, const double tau,
           const double dq, const double dphi, const double A)
{
    // Returns the D1 function that appears in the Scalar Singlet EOMs
    // Given q, phi, tau, dq and dphi as functions of A
    double kPhi = k(W0, xf, sc, ksc, kU1, kIR, k1, phi);
    double vf = Vf(W0, xf, sc, WIR, W1, phi, tau);
    double ans = -36864*exp(7*A) * kPhi * vf * pow(5 * q - 2 * dq, 4) * dphi / pow(q, 15);
    return ans;
}

double D2A(const double W0, const double xf, const double sc,
           const double WIR, const double W1, const double ksc,
           const double kU1, const double kIR, const double k1,
           const double q, const double phi, const double tau,
           const double dq, const double dtau, const double A)
{
    // Returns the D1 function that appears in the Scalar Singlet EOMs
    // Given q, phi, tau, dq and dphi as functions of A
    double kPhi = k(W0, xf, sc, ksc, kU1, kIR, k1, phi);
    double vf = Vf(W0, xf, sc, WIR, W1, phi, tau);
    double ans = -2304*exp(6*A) * pow(kPhi * vf,2) * pow(5 * q - 2 * dq, 2) * dtau / pow(q, 10);
    return ans;
}

double N1A(const double W0, const double xf, const double sc,
           const double ksc, const double kU1, const double kIR,
           const double VgIR, const double WIR, const double k1,
           const double W1, const double q, const double phi,
           const double tau, const double dq, const double dphi,
           const double dtau, const double A)
{
    // Returns the N1 function that appears in the Scalar Singlet EOMs
    // Given q, phi, tau, dq, dphi and dtau as functions of A
    double vg = Vg(sc, VgIR, phi);
    double dvg = dVgdPhi(sc, VgIR, phi);
    double vf = Vf(W0, xf, sc, WIR, W1, phi, tau);
    double dvfdphi = dVfdPhi(W0, xf, sc, WIR, W1, phi, tau);
    double dvfdtau = dVfdtau(W0,xf,sc,WIR,W1, phi, tau);
    double d2vfdphidtau = d2VfdPhidtau(W0, xf, sc, WIR, W1, phi, tau);
    double kPhi = k(W0, xf, sc, ksc, kU1, kIR, k1, phi);
    double dkdphi = dkdPhi(W0, xf, sc, ksc, kU1, kIR, k1, phi);
    double rA = RA(W0, xf, sc, ksc, kU1, kIR, k1, q, phi, dtau);
    double d1A = D1A(W0, xf, sc, WIR, W1, ksc, kU1, kIR, k1, q, phi, tau, dq, dphi, A);
    double par1 = q*q*(2*xf*vf+rA*vg)-rA*dphi*dphi;
    double par2 = -pow(1+rA*rA,2)*vf*dkdphi*dtau*dphi+4*pow(rA,4)*q*q*dvfdtau;
    double par3 = 4*vf*vf*(-3*(1+rA*rA)*dkdphi*pow(dtau*dphi,2)+rA*rA*pow(q,4)*(3*dvg-(4*rA*xf*vf+(1+rA*rA)*vg)*dphi))+12*pow(rA,3)*pow(q,4)*vf*(4*xf*vf+rA*vg)*dvfdphi-3*pow(rA,3)*pow(q,2)*dphi*(q*q*(2*xf*vf+rA*vg)-rA*dphi*dphi)*dvfdphi*dvfdphi;
    double par4 = 8*rA*rA*xf*pow(vf,3)*dkdphi*dtau-pow(rA,3)*pow(dphi,2)*dvfdtau*dvfdphi+2*(1+rA*rA)*vf*vf*dkdphi*dtau*(rA*vg-xf*dphi*dvfdphi)-rA*vf*dphi*(-8*rA*rA*dvfdtau+(1+rA*rA)*vg*dkdphi*dtau*dvfdphi+rA*rA*dphi*d2vfdphidtau);
    double par5 = rA*(1+rA*rA)*vf*dkdphi*dtau*pow(dphi,3)*dvfdphi+rA*rA*pow(q,4)*(2*xf*vf+rA*vg)*(dvfdphi*dvfdtau+vf*d2vfdphidtau) + q*q*par4;
    double ans = 16 * rA * pow(kPhi,4) * vf*vf * pow(dtau, 5) * pow(dphi, 3) + 3 * pow(q,4) * vf * dkdphi * par1 * par2;
    ans = ans - 96*pow(rA*kPhi,3)*q*q*vf*pow(dtau,3)*pow(dphi,2)*dvfdphi+4*rA*pow(kPhi*q,2)*dtau*par3;
    ans = ans + 12*rA*rA*kPhi*pow(q,4) * par5;
    ans = exp(7*A)*xf * dtau * ans / (d1A*pow(q,7)) ;
    return ans;
}

double N2A(const double W0, const double xf, const double sc,
           const double ksc, const double kU1, const double kIR,
           const double VgIR, const double WIR, const double k1,
           const double W1, const double q, const double phi,
           const double tau, const double dq, const double dphi,
           const double dtau, const double A)
{
    // Returns the N2 function that appears in the Scalar Singlet EOMs
    // Given q, phi, tau, dq, dphi and dtau as functions of A
    double vg = Vg(sc, VgIR, phi);
    double dvg = dVgdPhi(sc, VgIR, phi);
    double vf = Vf(W0, xf, sc, WIR, W1, phi, tau);
    double dvfdphi = dVfdPhi(W0, xf, sc, WIR, W1, phi, tau);
    double dvfdtau = dVfdtau(W0,xf,sc,WIR,W1, phi, tau);
    double d2vfdphidtau = d2VfdPhidtau(W0, xf, sc, WIR, W1, phi, tau);
    double d2vfdphi2 = d2Vfdphi2(W0, xf, sc, WIR, W1, phi, tau);
    double kPhi = k(W0, xf, sc, ksc, kU1, kIR, k1, phi);
    double dkdphi = dkdPhi(W0, xf, sc, ksc, kU1, kIR, k1, phi);
    double rA = RA(W0, xf, sc, ksc, kU1, kIR, k1, q, phi, dtau);
    double qA = QA(sc, VgIR, q, phi, dphi, A);
    double d2A = D2A(W0, xf, sc, WIR, W1, ksc, kU1, kIR, k1, q, phi, tau, dq, dtau, A);
    double par1 = -8*exp(2*A)*pow(rA*xf,2)*pow(q,4)*vf*vf*dvfdphi+2*exp(2*A)*rA*xf*vf*dkdphi*pow(dtau,2)*dphi*(vf*dphi+3*exp(4*A)*dvfdphi);
    par1 = par1 + pow(q,2)*(qA*vf*(vg*dkdphi*pow(dtau,2)+16*rA*rA*dphi)+4*qA*rA*rA*dphi*dphi*dvfdphi-2*exp(2*A)*xf*vf*(vf*((2*xf*vf+rA*vg)*dkdphi*dtau*dtau-4*(rA+pow(rA,3))*dphi)-2*rA*(2+rA*rA)*dphi*dphi*dvfdphi));
    double par2 = 6*qA*rA*rA*(1+rA*rA)*dkdphi*pow(dphi,2)+exp(2*A)*xf*vf*(3*(1+rA*rA)*(2*xf*vf+rA*vg)*dkdphi*pow(dtau,2)+12*(rA+pow(rA,3))*dkdphi*dphi*dphi-8*rA*pow(dphi,3))-12*exp(6*A)*rA*xf*pow(dphi,2)*dvfdphi;
    double par3 = -3*exp(2*A)*(rA+pow(rA,3))*xf*q*q*pow(vf,3)*dkdphi*dkdphi*pow(dtau,3)*pow(dphi,2)+pow(q,4)*vf*vf*dkdphi*dtau*par2;
    par3 = par3+pow(q,6)*(dkdphi*dtau*(3*(rA+pow(rA,3))*(qA*rA+2*exp(2*A)*xf*vf)*dvg+4*vf*vf*(qA*vg*dphi+4*exp(2*A)*xf*vf*(xf*vf+rA*vg)*dphi+3*exp(6*A)*pow(rA,4)*xf*(2*xf*vf+rA*vg)*dvfdphi))-4*rA*(qA*rA+2*exp(2*A)*xf*vf)*dphi*(dvfdtau*(2*vf*dphi-3*rA*rA*dvfdphi)+3*rA*rA*vf*d2vfdphidtau));
    double par4 = 6*exp(4*A)*qA*rA*rA*vf*dkdphi*pow(dphi*dtau,2)*dvfdphi+rA*pow(q,4)*vf*(-2*qA*rA*vf*dvg*dphi+qA*pow(rA,3)*(3*dvg+2*vg*dphi)*dvfdphi+2*exp(2*A)*xf*(-2*vf*vf*dvg*dphi+rA*rA*vf*(3*dvg+2*rA*rA*vg*dphi)*dvfdphi+3*pow(rA,3)*(2*xf*vf+rA*vg)*dvfdphi*dvfdphi));
    par4 = par4 +q*q*dphi*(-2*exp(2*A)*rA*xf*pow(vf,3)*(vg*dkdphi*pow(dtau,2)-16*dphi)-6*qA*pow(rA,4)*dphi*pow(dvfdphi,2)+6*pow(rA,3)*vf*dphi*(-exp(2*A)*(2+rA*rA)*xf*pow(dvfdphi,2)+qA*rA*d2vfdphi2)+pow(vf,2)*(qA*vg*dkdphi*pow(dtau,2)+4*exp(2*A)*pow(rA,2)*dphi*(4*vg+3*rA*xf*d2vfdphi2)));
    double ans = 12*rA*pow(q,4)*vf*(qA*rA+2*exp(2*A)*xf*vf)*dkdphi*dphi*(-vf*dkdphi*dtau*dphi+pow(q,2)*dvfdtau);
    ans = ans -2*pow(kPhi,3)*vf*pow(dtau,3)*dphi*par1 + kPhi*par3 + 2*pow(kPhi*q,2)*dtau*par4;
    ans = -exp(4*A)*ans/(d2A*pow(q,6));
    return ans;
}

struct eomScalarSingletFluctuations
{
    double W0, xf, sc, ksc, kU1, kIR;
    double VgIR, WIR, k1, W1;
    Spline_Interp<double> * qProfile, * phiProfile, * tauProfile;
    double msq;
    eomScalarSingletFluctuations(const double WW0, const double xxf, const double ssc,
                                 const double kksc, const double kkU1, const double kkIR,
                                 const double VVgIR, const double WWIR, const double kk1,
                                 const double WW1, Spline_Interp<double> q, Spline_Interp<double> phi,
                                 Spline_Interp<double> tau, const double m2)
                                 {
                                     W0 = WW0; xf = xxf; sc = ssc;
                                     ksc = kksc; kU1 = kkU1; kIR = kkIR;
                                     VgIR = VVgIR; WIR = WWIR; k1 = kk1;
                                     W1 = WW1; msq = m2;
                                     qProfile = &q;
                                     phiProfile = &phi;
                                     tauProfile = &tau;
                                 }
    void operator()(const state &X , state &dXdt , double A)
    {
        double q = qProfile->interp(A), phi = phiProfile->interp(A), tau = tauProfile->interp(A);
        double dq = qProfile->der1(A), dphi = phiProfile->der1(A), dtau = tauProfile->der1(A);
        double d2phi = phiProfile->der2(A);
        double n1a = N1A(W0, xf, sc, ksc, kU1, kIR, VgIR, WIR, k1, W1, q, phi, tau, dq, dphi, dtau, A);
        double n2a = N2A(W0, xf, sc, ksc, kU1, kIR, VgIR, WIR, k1, W1, q, phi, tau, dq, dphi, dtau, A);
        double kA = kSS(q, dq, dphi, d2phi);
        double pA = pSS(W0, xf, sc, ksc, kU1, kIR, WIR, k1, W1, q, phi, tau, dphi, dtau);
        double qA = qSS(W0, xf, sc, ksc, kU1, kIR, WIR, k1, W1, q, phi, tau, dphi, dtau);
        double tA = tSS(W0, xf, sc, ksc, kU1, kIR, k1, q, phi, dtau, A);
        double nA = nSS(W0, xf, sc, ksc, kU1, kIR, WIR, k1, W1, q, phi, tau, dq, dphi, dtau);
        // X = (zeta, xi, dzetadA, dxidA)
        dXdt[0] = X[2];
        dXdt[1] = X[3];
        dXdt[2] = n1a*(X[1]-X[0]) - exp(-2*A) * msq * q * q * X[0] - kA * X[2] - pA * X[3];
        dXdt[3] = n2a*(X[0]-X[1]) - msq * tA * X[1] - qA * X[2] - nA * X[3];  
    }
};

// Functions that define the IR boundary conditions
// zeta and xi are swaped relative to Jarvinen code
double XiIR(const double A, const double msq)
{
    return pow(sqrt(-A), 0.25*msq + 1);
}

double dXiIRdA(const double A, const double msq)
{
    return -0.5*(1+msq*0.25)*pow(-A,msq/8.0 - 0.5);
}

double XiIR2(const double A, const double msq)
{
    return pow(sqrt(-A), msq/6.0)/(1-msq/6.0);
}

double dXiIR2dA(const double A, const double msq)
{
    return -msq*pow(-A,-1+msq/12.)/(12*(1-msq/6.0));
}

double zetaIR(const double A, const double msq)
{
    return pow(sqrt(-A), msq/6.0);
}

double dzetaIRdA(const double A, const double msq)
{
    return -msq*pow(-A,msq/12. - 1) / 12.0;
}

vector<double> evolveScalarSingletFluctIRtoUV(const double sc, const double ksc, const double W0,
                                              const double kU1, const double kIR, const double VgIR,
                                              const double WIR, const double k1, const double W1,
                                              const double xf, const Spline_Interp<double> &qProf,
                                              const Spline_Interp<double> &phiProf,const Spline_Interp<double> &tauProf,
                                              const double m2, const string bc)
{
/* 
    Solves the Scalar Singlet EOMs from the IR to the UV.
    bc can either be 1 or 2 depending on the boundary conditions
    we are interested
*/
   // Setting up IR boundary conditions
   double AIR = -10, AUV = 4.0, h = 0.001;
    // Setting up the state object to start to solve the EOMs
    state X(4);
    // X = (zeta, xi, dzetadA, dxidA)
    if(bc == "1")
    {
        double zIR = zetaIR(AIR, m2);
        double dzetaIR = dzetaIRdA(AIR, m2);
        double xiIR = XiIR(AIR, m2) + XiIR2(AIR, m2);
        double dxiIR = dXiIRdA(AIR, m2) + dXiIR2dA(AIR, m2);
        X <<= zIR, xiIR, dzetaIR, dxiIR;
    }
    else if (bc == "2")
    {
        double zIR = -zetaIR(AIR, m2);
        double dzetaIR = -dzetaIRdA(AIR, m2);
        double xiIR = XiIR(AIR, m2) - XiIR2(AIR, m2);
        double dxiIR = dXiIRdA(AIR, m2) - dXiIR2dA(AIR, m2);
        X <<= zIR, xiIR, dzetaIR, dxiIR;
    }
    else
    {
        cout << "Boundary condition type input invalid" << endl;
        throw "error";
    }
    // Stepper definition
    dense_stepper stepper = make_dense_output( 1.0e-12 , 1.0e-12 , runge_kutta_dopri5< state >() );
    // Definition of the system of EOMs
    eomScalarSingletFluctuations eom(W0, xf, sc, ksc, kU1, kIR, VgIR, WIR, k1, W1, qProf, phiProf, tauProf, m2);
    integrate_const(stepper, eom, X , AIR , AUV, h);
    vector<double> sol(2);
    sol[0] = X[0]; sol[1] = X[1];
    // sol = zeta, xi
    return sol;
}

double detf(const double sc, const double ksc, const double W0,
            const double kU1, const double kIR, const double VgIR,
            const double WIR, const double k1, const double W1,
            const double xf, const Spline_Interp<double> &qProf,
            const Spline_Interp<double> &phiProf, const Spline_Interp<double> &tauProf,
            const double m2)
{
/* Computes the determinant of the matrix
   zeta1(AUV) xi1(AUV)
   zeta2(AUV) xi2(AUV)
   given by zeta1(AUV) * xi2(AUV) - xi1(AUV) * zeta2(AUV)
   Where the 1 and 2 correspond to solving the EOMs with bcs of type 1 and 2 respectively
*/
    // sol1 = zeta1, xi1
    vector<double> sol1 = evolveScalarSingletFluctIRtoUV(sc, ksc, W0, kU1, kIR, VgIR, WIR, k1, W1, xf, qProf, phiProf, tauProf, m2, "1");
    // sol2 = zeta2, xi2
    vector<double> sol2 = evolveScalarSingletFluctIRtoUV(sc, ksc, W0, kU1, kIR, VgIR, WIR, k1, W1, xf, qProf, phiProf, tauProf, m2, "2");
    // Return the determinant
    return sol1[0]*sol2[1] - sol1[1]*sol2[0];
}

vector<double> computeScalarSingletMasses(const double sc, const double ksc, const double W0,
                                          const double kU1, const double kIR, const double VgIR,
                                          const double WIR, const double k1, const double W1,
                                          const double xf, const Spline_Interp<double> &qProf,
                                          const Spline_Interp<double> &phiProf, const Spline_Interp<double> &tauProf,
                                          const int nmasses)
{
    // Computes the lowest nmasses of scalar singlet fluctuations
    // First we bracket the values
    // lbs and ubs will store the upper and lower bounds
    vector<double> lbs, ubs;
    double mmin=0.1, m = mmin, deltam = 0.03;
    double detf1 = detf(sc, ksc, W0, kU1, kIR, VgIR, WIR, k1, W1, xf, qProf, phiProf, tauProf, m*m);
    double detf2;
    int n = 0;
    while (n < nmasses)
    {
        detf2 = detf(sc, ksc, W0, kU1, kIR, VgIR, WIR, k1, W1, xf, qProf, phiProf, tauProf, pow(m+deltam,2));
        if (detf1*detf2 < 0)
        {
            lbs.push_back(m); ubs.push_back(m+deltam);
            m += deltam ;
            n += 1;
            detf1 = detf2;
        }
        else
        {
            m += deltam ;  
            detf1 = detf2;
        }
    }
    // Now that the mass bounds are obtained we can compute the masses
    vector<double> masses;
    function<double(double)> func = [&sc, &ksc, &W0, &kU1, &kIR, &VgIR, &WIR, &k1, &W1, &xf, &qProf, &phiProf, &tauProf] (double m) { return detf(sc, ksc, W0, kU1, kIR, VgIR, WIR, k1, W1, xf, qProf, phiProf, tauProf, m*m) ;} ;
    for(int i = 0; i < nmasses; i++)
    {
        double mass = zbrent(func, lbs[i], ubs[i], 1e-9, true);
        masses.push_back(mass);
    }
    return masses;
}

;
#endif