#include <iostream>
#include <fstream>
#include <functional>
#include <boost/numeric/ublas/assignment.hpp>
#include <boost/numeric/odeint.hpp>
#include "HolographicVQCD.h"
#include "schrodinger/common.h"
#include "methods/interpolation/Spline_Interp.hpp"

const std::vector<double> HVQCD::mrhos = {775.5, 1465, 1720, 1909, 2149, 2265};         // Non-Singlet Vector Mesons
const std::vector<double> HVQCD::ma1s = {1230, 1647, 1930, 2096, 2270};                 // Non-Singlet Axial Vector Mesons
const std::vector<double> HVQCD::mpis = {135, 1300, 1812, 2070, 2360};                  // Non-Singlet Pseudoscalar Mesons
const std::vector<double> HVQCD::a0s = {1474, 2025};                                    // Non-Singlet Scalar Mesons
const std::vector<double> HVQCD::mTG = {2150};                                          // Singlet Tensor glueballs
const std::vector<double> HVQCD::momegas = {782.65, 1420, 1670};                        // Singlet Vector Mesons
const std::vector<double> HVQCD::mf1s = {1281.9, 1426.4};                               // Singlet Axial Vector Mesons

// Ratios with the rho vector meson rho
const std::vector<double> HVQCD::RTG_rho = {2.7724};
const std::vector<double> HVQCD::Rrho_rho = {1.8891, 2.2179, 2.4616, 2.7711, 2.9207};
const std::vector<double> HVQCD::Ra1_rho = {1.5861, 2.1238, 2.48872, 2.70277, 2.92714};
const std::vector<double> HVQCD::Rpi_rho = {0.1741,1.6763,2.3366,2.6692, 3.0432};
const std::vector<double> HVQCD::Ra0_rho = {1.9007, 2.6112};
const std::vector<double> HVQCD::Romega_rho = {1.01, 1.83, 2.15};
const std::vector<double> HVQCD::Rf1_rho = {1.65, 1.84};


HVQCD::HVQCD(const double ssc, const double kksc, const double wwsc,
             const double WW0, const double ww0,
             const double kkU1, const double wwU1,
             const double VVgIR, const double WWIR, const double kkIR, const double wwIR,
             const double WW1, const double kk1, const double ww1,
             const double xxf, const double ttau0, 
             const double za, const double c,
             const bool add_m, const bool add_tg, const bool add_scals,
             const bool add_sing_v, const bool add_sing_av):
             Background(ssc, VVgIR), ksc(kksc), wsc(wwsc), W0(WW0), w0(ww0), 
             kU1(kkU1), wU1(wwU1), WIR(WWIR), kIR(kkIR), wIR(wwIR), W1(WW1),
             k1(kk1), w1(ww1), xf(xxf), tau0(ttau0),
             Za(za), ca(c),
             add_mass(add_m), add_tensor_glueball(add_tg), add_scalars(add_scals),
             add_singlet_vector(add_sing_v), add_singlet_axial(add_sing_av) {}

double HVQCD::Vf0l(const double l)
{
    // Definition of Vf0 as a function of lambda
    double coeff1 = (24 + 11 * W0 - 2 * W0 * xf) / (27 * M_PI * M_PI);
    double coeff2 = (24 * (857 - 46 * xf) + W0 * (4619 - 1714 * xf + 92 * xf * xf)) / (46656 * pow(M_PI, 4.0));
    double ans = W0 + coeff1 * l + coeff2 * l * l / (1 + sc * l / lambda0) + (3.0 / (16 * pow(M_PI,4.0))) * WIR * pow(l * sc, 2.0) * exp(- lambda0 / (sc * l)) * (1 + lambda0 * W1 / (sc * l)) ;
    return ans;
}

double HVQCD::dVf0dlambda(const double l)
{
    // Definition of dVf0/dlambda
    double v1 = (24 + 11 * W0 - 2 * W0 * xf)/(27*M_PI*M_PI);
    double v2 = (24 * (857 - 46 * xf) + W0 * (4619 - 1714 * xf + 92 * xf * xf))/(46656 * pow(M_PI,4) );
    double ans = v1 + v2 * l * lambda0 * (sc * l + 2 * lambda0) / pow(sc*l + lambda0,2);
    ans += 12 * exp(-lambda0/(sc*l))*WIR*(2*pow(sc*l,2)+sc*(1+W1)*l*lambda0+W1*pow(lambda0,2)) /(l*pow(lambda0,2));
    return ans;
}

double HVQCD::Vf0(const double phi)
{
    // Definition of Vf0 as a function of Phi, lambda = exp(Phi)
    double l = exp(phi);
    return Vf0l(l);
}

double HVQCD::dVf0dPhi(const double phi)
{
    // Definition of dVf0/dPhi
    double l = exp(phi);
    // First we compute dVf0/dlambda
    double ans = dVf0dlambda(l);
    // dVf0/dPhi = lambda * dVf0/dlambda
    return l * ans;
}

// Definition of the d2Vf0dPhi2 potential
double HVQCD::d2Vf0dPhi2(const double phi)
{
    // Definition of d2Vf0/dPhi2
    double coeff1 = (24 + 11 * W0 - 2 * W0 * xf) / (27 * M_PI * M_PI);
    double coeff2 = (24 * (857 - 46 * xf) + W0 * (4619 - 1714 * xf + 92 * xf * xf)) / (46656 * pow(M_PI, 4.0)); 
    double l = exp(phi);
    double denom = 1 + l * sc / lambda0;
    double exparg = lambda0/(sc*l);
    double ans = coeff1*l + 4*l*l*coeff2/denom + 2*pow(l,4)*sc*sc*coeff2/(pow(denom,3)*pow(lambda0,2))-5*pow(l,3)*sc*coeff2/(pow(denom,2)*lambda0);
    ans = ans -3*exp(-exparg+phi)*sc*W1*WIR*lambda0*(1+exparg)/(16*pow(M_PI,4));
    ans = ans -3*exp(-exparg+phi)*sc*W1*WIR*lambda0*(2+exparg)/(16*pow(M_PI,4))-3*exp(-exparg+phi)*sc*WIR*lambda0*(1+W1*exparg)/(16*pow(M_PI,4));
    ans = ans +3*exp(-exparg+2*phi)*sc*sc*WIR*(1+W1*exparg)*pow(2+exparg,2)/(16*pow(M_PI,4));
    return ans;
}

double HVQCD::d2Vf0dlambda2(const double l)
{
    // Definition of the d2Vf0/dlambda2 potential
    double ans = (d2Vf0dPhi2(log(l)) - dVf0dPhi(log(l)))/pow(l,2);
    return ans;
}

double HVQCD::Vfl(const double l, const double tau)
{
    // Definition of the Vf = Vf0 exp(-tau^2) potential
    return xf * Vf0l(l) * exp(-tau*tau);
}

double HVQCD::dVfldlambda(const double l, const double tau)
{
    // Definition of the dVfdl = dVf0dl exp(-tau^2) potential
    return xf * dVf0dlambda(l) * exp(-tau*tau);
}

double HVQCD::Vf(const double phi, const double tau)
{
    // Definition of the Vf = Vf0(Phi) exp(-tau^2) potential vs phi
    return xf * Vf0(phi) * exp(-tau * tau);
}

double HVQCD::dVfdPhi(const double phi, const double tau)
{
    // Definition of the dVfdPhi = dVf0dPhi exp(-tau^2) potential
    return xf * dVf0dPhi(phi) * exp(-tau*tau);
}

// Declaration of the dVfdtau = Vf0 (-2tau)exp(-tau^2) potential
double HVQCD::dVfdtau(const double phi, const double tau)
{
    // Definition of the dVfdtau = Vf0 (-2tau)exp(-tau^2) potential
    return xf * Vf0(phi) * (-2 * tau) * exp(-tau*tau);
}

double HVQCD::d2VfdPhi2(const double phi, const double tau)
{
    // Definition of d2Vf/dPhi2 potential
    return xf * d2Vf0dPhi2(phi) * exp(-tau*tau);
}

double HVQCD::d2VfdPhidtau(const double phi, const double tau)
{
    // Definition of d2Vf/dPhidtau = dVf0dPhi (-2tau)exp(-tau^2)
    return xf * dVf0dPhi(phi) * (-2 * tau) * exp(-tau*tau);
}

double HVQCD::klambda(const double l)
{
    // Definition of the k potential as a function of lambda = exp(Phi)
    double ans = 1.5 - W0 * xf / 8;
    ans = ans * (1 + sc * kU1 * l / lambda0 + kIR * exp(-lambda0/(ksc * l)) * (1 + lambda0 * k1 /(ksc * l)) * pow(ksc * l / lambda0, 4.0/3) / sqrt(log(1 + ksc * l / lambda0 ))) ;
    return 1 / ans ;
}

double HVQCD::k(const double phi)
{
    // Definition of the k potential as a function of Phi
    double l = exp(phi);
    return klambda(l);
}

// Definition of dk/dl potential
double HVQCD::dkdlambda(const double l)
{
    // Computes dk/dlambda
    double k0 = 1.5 - W0 * xf / 8.0;
    double kfactor = k0 * pow( 1 + kU1 * sc * l / lambda0 + exp(-lambda0/(ksc * l)) * kIR * pow(ksc * l, 4.0/3) * (1+ k1 * lambda0 /(ksc * l))/(16 * pow(M_PI, 8.0/3) * sqrt(log(1+ksc*l/lambda0))), 2.0) ;
    double ans = kU1 * sc / lambda0 - exp(-lambda0/(ksc*l)) * kIR * ksc * pow(ksc * l, 4.0/3) * (1+ k1 * lambda0 /(ksc * l)) / (32 * pow(M_PI, 8.0/3) * (1+ksc * l / lambda0) * lambda0 * pow(log(1+ksc*l/lambda0),1.5));
    ans -= exp(-lambda0/(ksc*l)) * k1 * kIR * pow(ksc * l, 4.0/3) * lambda0 / (16 * ksc * pow(M_PI,8.0/3) * l * l * sqrt(log(1+ksc*l/lambda0)));
    ans += exp(-lambda0/(ksc*l)) * kIR * ksc * pow(ksc * l, 1.0/3) * (1+k1*lambda0/(ksc * l)) / (12 * pow(M_PI,8.0/3) * sqrt(log(1+ksc*l/lambda0)));
    ans += exp(-lambda0/(ksc*l)) * kIR * pow(ksc * l, 4.0/3) * lambda0 * (1+k1*lambda0/(ksc * l)) / (16 * ksc * pow(M_PI, 8.0/3) * l * l * sqrt(log(1+ ksc*l / lambda0)));
    ans = - ans / kfactor ;
    return ans;
}

double HVQCD::dkdPhi(const double phi)
{
    // Computes dk/dPhi = lambda dk/dlambda
    double l = exp(phi);
    // Compute dk/dlambda
    double ans = dkdlambda(l);
    // Return lambda * dk/dlambda = dk/dPhi
    return l * ans;
}

double HVQCD::d2kdPhi2(const double phi)
{
    /* Definition of d2k/dPhi2 potential
       Given that k(Phi) = 1 / f(Phi), d2k/dPhi2 = 2 dkdPhi^2/k - k^2 f''(Phi)
    */
    double kPhi = k(phi);
    double dkPhi = dkdPhi(phi);
    double l = exp(phi);
    double arg = ksc*l/lambda0;
    double logarg = log(1+arg);
    double d2f = kU1*sc*l/lambda0 + 0.75*exp(-1/arg)*kIR*pow(ksc*l,2)*pow(arg,4./3)*(1+k1/arg)/(pow((1+arg)*lambda0,2)*pow(logarg,2.5)) + exp(-1/arg)*k1*kIR*pow(arg,4./3)/((1+arg)*pow(logarg,1.5));
    d2f += - exp(-1/arg)*kIR*pow(arg,4./3)*(1+k1/arg)/(2*(1+arg)*pow(logarg,1.5)) - (4./3)*exp(-1/arg)*kIR*pow(ksc*l,2)*pow(arg,1./3)*(1+k1/arg)/((1+arg)*pow(lambda0,2)*pow(logarg,1.5));
    d2f += exp(-1/arg)*kIR*pow(ksc*l,2)*pow(arg,4./3)*(1+k1/arg)/(2*pow((1+arg)*lambda0,2)*pow(logarg,1.5))-exp(-1/arg)*kIR*ksc*l*pow(arg,4./3)*(1+1/arg)*(1+k1/arg)/(2*(1+arg)*lambda0*pow(logarg,1.5));
    d2f += -(8./3)*exp(-1/arg)*k1*kIR*pow(arg,1./3)/sqrt(logarg)-exp(-1/arg)*k1*kIR*pow(arg,4./3)*pow(lambda0/(ksc*l),2)/sqrt(logarg)-exp(-1/arg)*k1*kIR*pow(arg,4./3)*lambda0*(-1+1/arg)/(ksc*l*sqrt(logarg));
    d2f += (4./3)*exp(-1/arg)*kIR*pow(arg,1./3)*(1+k1/arg)/sqrt(logarg) +(4./9)*exp(-1/arg)*kIR*pow(ksc*l/lambda0,2)*(1+k1/arg)/(pow(arg,2./3)*sqrt(logarg));
    d2f += exp(-1/arg)*kIR*pow(arg,4./3)*lambda0*(-1+1/arg)*(1+k1/arg)/(ksc*l*sqrt(logarg)) +(4./3)*exp(-1/arg)*kIR*ksc*l*pow(arg,1./3)*(1+1/arg)*(1+k1/arg)/(lambda0*sqrt(logarg));
    d2f = (1.5 - W0*xf/8)*d2f;
    double ans = 2 * pow(dkPhi,2)/kPhi - pow(kPhi,2)*d2f; 
    return ans;
}

double HVQCD::d2kdlambda2(const double l)
{
    // Compute d2k/dlambda2
    double ans = (d2kdPhi2(log(l)) - dkdPhi(log(l)))/pow(l,2);
    return ans;
}

double HVQCD::w(const double phi)
{
    // Returns w(Phi) potential
    double l = exp(phi);
    double ans = w0 * (1 + sc * wU1 * l / (lambda0 * (1 + sc * l / lambda0)) + wIR * exp(-lambda0 /(wsc * l)) * (1+lambda0*w1/(wsc * l)) * pow(wsc * l / lambda0, 4.0/3) / log(1 + wsc * l / lambda0));
    return 1 / ans;
}

double HVQCD::dwdPhi(const double phi)
{
    /* Definition of the dw/dPhi potential
       Given that w(Phi) = 1 / f(Phi), dw/dPhi = -f'(Phi)w(Phi)^2
    */
    double wPhi = w(phi);
    double l = exp(phi);
    double arg = wsc*l/lambda0;
    double logarg = log(1+arg);
    double df = -pow(sc*l/lambda0,2)*wU1/pow(1+sc*l/lambda0,2)+sc*wU1*l/((1+sc*l/lambda0)*lambda0)-exp(-1/arg)*wIR*wsc*l*pow(arg,4./3)*(1+w1/arg)/((1+arg)*lambda0*pow(logarg,2));
    df += -exp(-1/arg)*w1*wIR*pow(arg,4./3)*lambda0/(wsc*l*logarg)+(4./3)*exp(-1/arg)*wIR*wsc*l*pow(arg,1./3)*(1+w1/arg)/(lambda0*logarg)+exp(-1/arg)*wIR*pow(arg,4./3)*lambda0*(1+w1/arg)/(wsc*l*logarg);
    df = w0 * df;
    double ans = - df*pow(wPhi,2);
    return ans;
}

double HVQCD::dwdlambda(const double l)
{
    // Returns dw/dlambda = dw/dPhi / lambda
    double ans = dwdPhi(log(l));
    return ans/l;
}

// Definition of the d2w/dPhi2 potential
double HVQCD::d2wdPhi2(const double phi)
{
    // Given that w(Phi) = 1 / f(Phi), d2w/dPhi2 = 2 dwdPhi^2/w - w^2 f''(Phi)
    double wPhi = w(phi);
    double dwPhi = dwdPhi(phi);
    double l = exp(phi);
    double arg = wsc*l/lambda0;
    double logarg = log(1+arg);
    double d2f = 2*pow(sc*l/(lambda0*(1+sc*l/lambda0)),3)*wU1-3*pow(sc*l/((1+sc*l/lambda0)*lambda0),2)*wU1+sc*wU1*l/((1+sc*l/lambda0)*lambda0)+2*exp(-1/arg)*wIR*pow(wsc*l/lambda0,2)*pow(arg,4./3)*(1+w1/arg)/(pow(1+arg,2)*pow(logarg,3));
    d2f += 2*exp(-1/arg)*w1*wIR*pow(arg,4./3)/((1+arg)*pow(logarg,2))-exp(-1/arg)*wIR*pow(arg,4./3)*(1+w1/arg)/((1+arg)*pow(logarg,2))-(8./3)*exp(-1/arg)*wIR*pow(wsc*l/lambda0,2)*pow(arg,1./3)*(1+w1/arg)/((1+arg)*pow(logarg,2));
    d2f += exp(-1/arg)*wIR*pow(wsc*l/lambda0,2)*pow(arg,4./3)*(1+w1/arg)/(pow(1+wsc*l/lambda0,2)*pow(logarg,2)) -exp(-1/arg)*wIR*wsc*l*pow(arg,4./3)*(1+1/arg)*(1+w1/arg)/((1+arg)*lambda0*pow(logarg,2));
    d2f += -(8./3)*exp(-1/arg)*w1*wIR*pow(arg,1./3)/logarg-exp(-1/arg)*w1*wIR*pow(arg,4./3)*pow(lambda0/(wsc*l),2)/logarg-exp(-1/arg)*w1*wIR*pow(arg,4./3)*lambda0*(-1+1/arg)/(wsc*l*logarg);
    d2f += (4./3)*exp(-1/arg)*wIR*pow(arg,1./3)*(1+w1/arg)/logarg+(4./9)*exp(-1/arg)*wIR*pow(wsc*l/lambda0,2)*(1+w1/arg)/(pow(arg,2./3)*logarg);
    d2f += exp(-1/arg)*wIR*pow(arg,4./3)*lambda0*(-1+1/arg)*(1+w1/arg)/(wsc*l*logarg)+(4./3)*exp(-1/arg)*wIR*wsc*l*pow(arg,1./3)*(1+1/arg)*(1+w1/arg)/(lambda0*logarg);
    d2f = w0 * d2f;
    double ans = 2*pow(dwPhi,2)/wPhi -pow(wPhi,2)*d2f;
    return ans;
}

double HVQCD::G(const double q, const double phi, const double dt)
{
    return sqrt(1 + k(phi) * pow(dt / q, 2.0));
}

double HVQCD::dG(const double q, const double phi, const double dq, const double dphi, const double dt, const double d2t)
{
    // Returns dG/dA
    double kPhi = k(phi);
    double dkPhi = dkdPhi(phi);
    double g = G(q, phi, dt);
    double ans = (-kPhi*dq*pow(dt,2)/q+dkPhi*pow(dt,2)*dphi/2+kPhi*dt*d2t)/(q*q*g);
    return ans;
}

double HVQCD::d2G(const double q, const double phi, const double dq,
                  const double dphi, const double dtau, const double d2q,
                  const double d2phi, const double d2tau, const double d3tau)
{
    // Returns d2G/dA2
    double g = G(q, phi, dtau);
    double kPhi = k(phi);
    double dkPhi = dkdPhi(phi);
    double d2kPhi = d2kdPhi2(phi);
    double ans = - pow(kPhi*dq*dtau*dtau,2)/pow(q*q*g,3)+3*kPhi*pow(dq*dtau/(q*q),2)/g+kPhi*dkPhi*dq*pow(dtau,4)*dphi/(pow(q,5)*pow(g,3))-2*dkPhi*dq*pow(dtau,2)*dphi/(pow(q,3)*g);
    ans += -pow(dkPhi*dtau*dtau*dphi/(q*q),2)/(4*pow(g,3))+pow(dtau*dphi/q,2)*d2kPhi/(2*g)-kPhi*pow(dtau,2)*d2q/(pow(q,3)*g)+2*pow(kPhi,2)*dq*pow(dtau,3)*d2tau/(pow(q,5)*pow(g,3));
    ans += -4*kPhi*dq*dtau*d2tau/(pow(q,3)*g)-kPhi*dkPhi*pow(dtau,3)*dphi*d2tau/(pow(q,4)*pow(g,3))+2*dkPhi*dtau*dphi*d2tau/(pow(q,2)*g);
    ans += -pow(kPhi*dtau*d2tau/(q*q),2)/pow(g,3)+kPhi*pow(d2tau/q,2)/g+dkPhi*pow(dtau/q,2)*d2phi/(2*g)+kPhi*dtau*d3tau/(q*q*g);
    return ans;
}

double HVQCD::Z(const double l)
{
    // Returns Z(lambda) = Za(1 + ca * l ^ 4)
    return Za*(1 + ca * std::pow(l,4));
}

double HVQCD::dudA(const double q, const double phi,
                   const double dtau, const double A)
{
    // Returns du/dA = G(A) q(A) exp(-A)
    return G(q, phi, dtau)*q*exp(-A);
}

double HVQCD::dtauYangMills1(const double q, const double phi, const double tau, const double dphi)
{
    /* 
       Returns the first derivative of tau in the YM regime
       Valid for tau > tcut
    */
    double vf0 = Vf0(phi);
    double dvf0 = dVf0dPhi(phi);
    double kphi = k(phi);
    double dk = dkdPhi(phi);
    double ans = - 4 * q * q * vf0 * tau;
    ans = ans / (8*vf0 * kphi + 2 * kphi * dvf0 * dphi + vf0 * dk * dphi) ;
    return ans;
}

double HVQCD::d2tauYangMills1(const double q, const double phi, const double tau, const double dq, const double dphi, const double d2phi)
{
    /* 
       Returns the second derivative of tau in the YM regime
       Valid for tau > tcut
    */
    double vf0 = Vf0(phi), dvf0 = dVf0dPhi(phi);
    double d2vf0 = d2Vf0dPhi2(phi);
    double kPhi = k(phi), dkPhi = dkdPhi(phi);
    double d2kPhi = d2kdPhi2(phi);
    double num = 16*pow(q,4)*pow(vf0,2)*tau-64*q*pow(vf0,2)*kPhi*tau*dq-16*q*vf0*kPhi*tau*dq*dvf0*dphi;
    num += 32*pow(q*vf0,2)*tau*dkPhi*dphi-8*q*pow(vf0,2)*tau*dq*dkPhi*dphi-8*pow(q,2)*kPhi*tau*pow(dvf0*dphi,2);
    num += 8*pow(q,2)*vf0*tau*dvf0*dkPhi*pow(dphi,2)+8*q*q*vf0*kPhi*tau*pow(dphi,2)*d2vf0+4*pow(q*vf0*dphi,2)*tau*d2kPhi;
    num += 8*pow(q,2)*vf0*kPhi*tau*dvf0*d2phi+4*pow(q*vf0,2)*tau*dkPhi*d2phi;
    double denom = pow(8*vf0*kPhi+(2*kPhi*dvf0+vf0*dkPhi)*dphi,2);
    double d2tau = num/denom;
    return d2tau;
}

double HVQCD::d2tauYM2dA2(const double q, const double phi, const double tau, const double dq, const double dphi, const double dtau)
{
    double vf0 = Vf0(phi);
    double dvf0 = dVf0dPhi(phi);
    double kl = k(phi);
    double dkl = dkdPhi(phi);
    double ans = - 2 * q * q * tau / kl - 4 * dtau + dq * dtau / q - 2 * tau * pow(dtau,2) - 4 * kl * pow(dtau,3)/pow(q,2) - dvf0 * dtau * dphi/vf0;
    ans = ans - dkl * dtau * dphi / kl - kl * dvf0 * pow(dtau,3) * dphi / (pow(q,2)*vf0) - dkl * pow(dtau,3) * dphi / (2*q*q);
    return ans;
}

double HVQCD::d3tauYM(const double q, const double phi, const double tau,
                      const double dq, const double dphi, const double dtau,
                      const double d2q, const double d2phi, const double d2tau)
{
    // Returns d3tau/dA3 in the Yang-Mills regime
    double kPhi = k(phi);
    double dkPhi = dkdPhi(phi);
    double d2kPhi = d2kdPhi2(phi);
    double vf0 = Vf0(phi);
    double dvf0 = dVf0dPhi(phi);
    double d2vf0 = d2Vf0dPhi2(phi);
    double ans = -4*q*tau*dq/kPhi-2*q*q*dtau/kPhi-dq*dq*dtau/(q*q)-2*pow(dtau,3)+8*kPhi*dq*pow(dtau,3)/pow(q,3)+2*q*q*tau*dkPhi*dphi/pow(kPhi,2);
    ans += -4*dkPhi*pow(dtau,3)*dphi/pow(q,2)+dq*dkPhi*pow(dtau/q,3)*dphi+2*kPhi*dq*dvf0*pow(dtau/q,3)*dphi/vf0+pow(dkPhi*dphi/kPhi,2)*dtau;
    ans += pow(dvf0*dphi/vf0,2)*dtau-dvf0*dkPhi*pow(dtau,3)*pow(dphi/q,2)/vf0+kPhi*pow(dvf0*dphi/(q*vf0),2)*pow(dtau,3)-dtau*pow(dphi,2)*d2kPhi/kPhi;
    ans += -pow(dtau,3)*pow(dphi/q,2)*d2kPhi/2. + dtau*d2q/q - dtau*dphi*dphi*d2vf0/vf0-kPhi*pow(dtau,3)*pow(dphi/q,2)*d2vf0/vf0;
    ans += -4*d2tau + dq*d2tau/q -4*tau*dtau*d2tau -12*kPhi*pow(dtau/q,2)*d2tau-dkPhi*dphi*d2tau/kPhi;
    ans += -dvf0*dphi*d2tau/vf0 -3*dkPhi*pow(dtau/q,2)*dphi*d2tau/2. - 3*kPhi*dvf0*pow(dtau/q,2)*dphi*d2tau/vf0;
    ans += -dkPhi*dtau*d2phi/kPhi - dvf0*dtau*d2phi/vf0 - dkPhi*pow(dtau,3)*d2phi/(2.*q*q) - kPhi*dvf0*pow(dtau,3)*d2phi/(q*q*vf0);
    return ans;
}

void HVQCD::eomYangMills1(const state &X, state &dXdA, const double A)
{
    // X = q, Phi, tau, z, u
    dXdA[0] = dqYM(X[0], X[1]);
    dXdA[1] = dPhiYM(X[0], X[1]);
    dXdA[2] = dtauYangMills1(X[0], X[1], X[2], dXdA[1]);
    dXdA[3] = dzdA(X[0], A);
    dXdA[4] = dudA(X[0], X[1], dXdA[2], A);
}

void HVQCD::eomYangMills2(const state &X, state &dXdA, const double A)
{
    // X = (q, Phi, tau, dtau, z, u)
    dXdA[0] = dqYM(X[0], X[1]);
    dXdA[1] = dPhiYM(X[0], X[1]);
    dXdA[2] = X[3];
    dXdA[3] = d2tauYM2dA2(X[0], X[1], X[2], dXdA[0], dXdA[1], X[3]); 
    dXdA[4] = dzdA(X[0], A);
    dXdA[5] = dudA(X[0], X[1], X[3], A);
}

void HVQCD::observerYangMills2(const state &X , const double A)
{
    // X = (q, phi, tau, dtau/dA, z, u)
    double dq = dqYM(X[0], X[1]);
    double dphi = dPhiYM(X[0], X[1]);
    double d2q = d2qYM(X[0], X[1]);
    double d2phi = d2PhiYM(X[0], X[1]);
    double d2tau = d2tauYM2dA2(X[0], X[1], X[2], dq, dphi, X[3]);
    qYM2.push_back(X[0]);
    PhiYM2.push_back(X[1]);
    tauYM2.push_back(X[2]);
    dqYM2.push_back(dq);
    dPhiYM2.push_back(dphi);
    dtauYM2.push_back(X[3]);
    d2qYM2.push_back(d2q);
    d2PhiYM2.push_back(d2phi);
    d2tauYM2.push_back(d2tau);
    //d3taus.push_back(d3tauYM(W0, xf, sc, ksc, kU1, WIR, kIR, k1, W1, X[0], X[1], X[2], dq, dphi, X[3], d2q, d2phi, d2tau));
    d3tauYM2.push_back(d3tauCoupled(X[0], X[1], X[2], dq, dphi, X[3], d2phi, d2tau));
    AYM2.push_back(A);
    zYM2.push_back(X[4]);
    uYM2.push_back(X[5]);
}

double HVQCD::d3tauCoupled(const double q, const double phi, const double tau, const double dq, const double dphi,
                           const double dtau, const double d2phi, const double d2tau)
{
    // Returns d3tau/dA3 in the coupled regime
    double kPhi = k(phi);
    double dkPhi = dkdPhi(phi);
    double d2kPhi = d2kdPhi2(phi);
    double vf0 = Vf0(phi);
    double dvf0 = dVf0dPhi(phi);
    double d2vf0 = d2Vf0dPhi2(phi);
    double vf = Vf(phi, tau);
    double dvfdphi = dVfdPhi(phi, tau);
    double dvfdtau = dVfdtau(phi, tau);
    double g = G(q, phi, dtau);
    double dg = dG(q, phi, dq, dphi, dtau, d2tau);
    double ans = -4*q*tau*dq/kPhi-2*q*q*dtau/kPhi-2*pow(dtau,3) -kPhi*vf*dg*pow(dtau,3)/(6*pow(g,2))+8*kPhi*dq*pow(dtau/q,3);
    ans += 2*pow(q/kPhi,2)*tau*dkPhi*dphi -4*dkPhi*pow(dtau,3)*dphi/pow(q,2) +vf*dkPhi*pow(dtau,3)*dphi/(6*g);
    ans += dq*dkPhi*pow(dtau/q,3)*dphi+2*kPhi*dq*dvf0*pow(dtau/q,3)*dphi/vf0+pow(dkPhi*dphi/kPhi,2)*dtau;
    ans += pow(dvf0*dphi/vf0,2)*dtau-dvf0*dkPhi*pow(dtau,3)*pow(dphi/q,2)/vf0+kPhi*pow(dvf0*dphi/(q*vf0),2)*pow(dtau,3);
    ans += -dtau*pow(dphi,2)*d2kPhi/kPhi-pow(dtau,3)*pow(dphi/q,2)*d2kPhi/2 -dtau*pow(dphi,2)*d2vf0/vf0-kPhi*pow(dtau,3)*pow(dphi/q,2)*d2vf0/vf0;
    ans += -4*d2tau-4*tau*dtau*d2tau-12*kPhi*pow(dtau/q,2)*d2tau+vf*kPhi*pow(dtau,2)*d2tau/(2*g)-dkPhi*dphi*d2tau/kPhi;
    ans += -dvf0*dphi*d2tau/vf0-3*dkPhi*pow(dtau/q,2)*dphi*d2tau/2-3*kPhi*dvf0*pow(dtau/q,2)*dphi*d2tau/vf0+(4./9)*pow(dphi,2)*d2tau;
    ans += -dkPhi*dtau*d2phi/kPhi-dvf0*dtau*d2phi/vf0 -dkPhi*pow(dtau,3)*d2phi/(2*q*q)-kPhi*dvf0*pow(dtau,3)*d2phi/(q*q*vf0);
    ans += (8./9)*dtau*dphi*d2phi+kPhi*pow(dtau,4)*dvfdtau/(6*g)+kPhi*pow(dtau,3)*dphi*dvfdphi/(6*g);
    return ans;
}

double HVQCD::dqCoupled(const double q, const double phi, const double tau, const double dphi, const double dtau)
{
    // Returns dq/dA form the coupled EOMs of q, Phi and tau
   double kPhi = k(phi);
   double ans = (4.0 / 9) * q * pow(dphi,2) + q * Vf(phi, tau) * kPhi * pow(dtau,2) / (6.0 * sqrt(1 + kPhi * pow(dtau/q,2.0)));
   return ans;
}

double HVQCD::d2qCoupled(const double q, const double phi, const double tau,  const double dq, const double dphi,
                         const double dtau, const double d2tau)
{
    // Returns d2q/dA2 form the coupled EOMs of q, Phi and tau
   double vg = Vg(phi);
   double dvg = dVgdPhi(phi);
   double vf = Vf(phi, tau);
   double dvfdphi = dVfdPhi(phi, tau);
   double dvfdtau = dVfdtau(phi, tau);
   double kPhi = k(phi);
   double dkPhi = dkdPhi(phi);
   double g = G(q, phi, dtau);
   double dg = dG(q, phi, dq, dphi, dtau, d2tau);
   double d2q = -pow(q,3)*vf*dg/(3*g*g)+4*dq + q*q*vf*dq/g-q*q*vg*dq-kPhi*q*vf*dg*pow(dtau,2)/(6*g*g);
   d2q += kPhi*vf*dq*pow(dtau,2)/(6*g) - pow(q,3)*dvg*dphi/3 + q*vf*dkPhi*pow(dtau,2)*dphi/(6*g);
   d2q += kPhi*q*vf*dtau*d2tau/(3*g) + pow(q,3)*dtau*dvfdtau/(3*g)+kPhi*q*pow(dtau,3)*dvfdtau/(6*g);
   d2q += pow(q,3)*dphi*dvfdphi/(3*g)+kPhi*q*pow(dtau,2)*dphi*dvfdphi/(6*g);
   return d2q;
}

double HVQCD::d2PhiCoupled(const double q, const double phi, const double tau, const double dphi, const double dtau)
{
    // Returns d2Phi/dA2 from the coupled EOMs of q, lambda and tau
    double vg = Vg(phi);
    double dvg = dVgdPhi(phi);
    double vf = Vf(phi, tau);
    double dvfdphi = dVfdPhi(phi, tau);
    double kl = k(phi);
    double dkdphi = dkdPhi(phi);
    double g = sqrt(1 + kl * pow(dtau/q,2));
    double ans = -(3./8) * q * q * dvg + 9./dphi + (3./4) * q * q * (vf/g - vg)/dphi;
    ans = ans - 5 * dphi + vf * kl * pow(dtau,2) * dphi / (6*g) + (4./9) * pow(dphi,3);
    ans = ans + (3./8)*q*q*dvfdphi/g + (3./8)*kl*dtau*dtau*dvfdphi/g;
    ans = ans + (3./16)*vf*dtau*dtau*dkdphi/g;
    return ans;
}

double HVQCD::dPhiConstr(const double q, const double phi, const double tau, const double dtau)
{
    double vg = Vg(phi);
    double vf = Vf(phi, tau);
    double kPhi = k(phi);
    double ans = -0.5 * sqrt(3) * sqrt(12 -pow(q,2) * vg + pow(q,2) * vf / sqrt(1 + kPhi * dtau * dtau / pow(q,2)));
    return ans;
}

double HVQCD::d2tauCoupled(const double q, const double phi, const double tau, const double dphi, const double dtau)
{
    /*
        Returns d2tau from the coupled EOMs of q, lambda and tau
    */
    // X = (q, phi, tau, dphi, dtau)
    double dvf0 = dVf0dPhi(phi);
    double vf0 = Vf0(phi);
    double vf = Vf(phi, tau);
    double kl = k(phi);
    double dkl = dkdPhi(phi);
    double g = sqrt(1 + kl * pow(dtau/q,2.0));
    double ans = - 2 * q * q * tau / kl - 4 * dtau - 2 * tau * pow(dtau,2.0) - 4 * kl * pow(dtau,3.0) / pow(q,2.0) + kl * vf * pow(dtau,3.0) / (6.0 * g); 
    ans = ans - dvf0 * dphi * dtau / vf0 - kl * dvf0 * dphi * pow(dtau,3.0) / (pow(q,2.0) * vf0) + (4./9) * pow(dphi,2.0) * dtau;
    ans = ans  - dkl * dphi * dtau / kl - 0.5 * dkl * dphi * pow(dtau, 3.0) / pow(q,2.0) ;
    return ans;
}

void HVQCD::jacobian(const state &X , matrix_type &jac , const double A, state &dfdt)
{
    /* Computes the jacobian of the system of the coupled EOMs
       X = (q, phi, tau, dphi, dtau, z) = (X1, X2, X3, X4, X5, X6)
    */
    double kPhi = k(X[1]);
    double dkPhi = dkdPhi(X[1]);
    double d2kPhi = d2kdPhi2(X[1]);
    double vg = Vg(X[1]);
    double dvg = dVgdPhi(X[1]);
    double d2vg = d2VgdPhi2(X[1]);
    double vf0 = Vf0(X[1]);
    double dvf0 = dVf0dPhi(X[1]);
    double d2vf0 = d2Vf0dPhi2(X[1]);
    double vf = Vf(X[1], X[2]);
    double dvfdphi = dVfdPhi(X[1], X[2]);
    double d2vfdphi2 = d2VfdPhi2(X[1], X[2]);
    double dvfdtau = dVfdtau(X[1], X[2]);
    double d2vfdphidtau = d2VfdPhidtau(X[1], X[2]);
    double g = sqrt(1 + kPhi * pow(X[4]/X[0],2.0));
    jac(0,0) = (4./9) * pow(X[3],2) + pow(X[4],4)*pow(kPhi/X[0],2) * vf / (6*pow(g,3)) + pow(X[4],2)*kPhi*vf/(6*g);
    jac(0,1) = -pow(X[4],4)*kPhi*vf*dkPhi/(12*X[0]*pow(g,3)) + X[0]*pow(X[4],2)*vf*dkPhi/(6*g) + X[0]*pow(X[4],2)*kPhi*dvfdphi/(6*g);
    jac(0,2) = X[0]*pow(X[4],2)*kPhi*dvfdtau/(6*g);
    jac(0,3) = (8./9)*X[0]*X[3];
    jac(0,4) = -pow(X[4],3)*pow(kPhi,2)*vf/(6*X[0]*pow(g,3))+X[0]*X[4]*kPhi*vf/(3*g);
    jac(0,5) = 0; jac(0,6) = 0;
    jac(1,0) = 0; jac(1,1) = 0; jac(1, 2) = 0; jac(1, 3) = 1; jac(1,4) = 0; jac(1,5) = 0; jac(1,6) = 0;
    jac(2,0) = 0; jac(2,1) = 0; jac(2, 2) = 0; jac(2, 3) = 0; jac(2,4) = 1; jac(2,5) = 0; jac(2,6) = 0;
    jac(3,0) = 0.75*pow(X[4],2)*kPhi*vf/(X[0]*X[3]*pow(g,3))+X[3]*pow(X[4],4)*pow(kPhi,2)*vf/(6*pow(X[0]*g,3));
    jac(3,0) += 1.5*X[0]*vf/(X[3]*g)-1.5*X[0]*vg/X[3]+(3./16)*pow(X[4],4)*kPhi*vf*dkPhi/pow(X[0]*g,3)-0.75*X[0]*dvg;
    jac(3,0) += (3./8)*pow(X[4],2)*kPhi*dvfdphi/(X[0]*pow(g,3))+(3./8)*pow(X[4],4)*pow(kPhi,2)*dvfdphi/(pow(X[0]*g,3)) + 0.75*X[0]*dvfdphi/g;
    jac(3,1) = -(3./8)*pow(X[4],2)*vf*dkPhi/(X[3]*pow(g,3)) -(1./12)*X[3]*pow(X[4],4)*kPhi*vf*dkPhi/(pow(X[0],2)*pow(g,3)) + X[3]*pow(X[4],2)*vf*dkPhi/(6*g) -(3./32)*pow(X[4],4)*vf*pow(dkPhi/X[0],2)/pow(g,3);
    jac(3,1) += -0.75*X[0]*X[0]*dvg/X[3] + (3./16)*pow(X[4],2)*vf*d2kPhi/g -(3./8)*pow(X[0],2)*d2vg;
    jac(3,1) += 0.75*pow(X[0],2)*dvfdphi/(X[3]*g) + X[3]*pow(X[4],2)*kPhi*dvfdphi/(6*g) -(3./16)*pow(X[4],2)*dkPhi*dvfdphi/pow(g,3) - (3/16.)*pow(X[4],4)*kPhi*dkPhi*dvfdphi/(pow(X[0],2)*pow(g,3));
    jac(3,1) += (9./16)*pow(X[4],2)*dkPhi*dvfdphi/g + (3./8)*pow(X[0],2)*d2vfdphi2/g + (3./8)*pow(X[4],2)*kPhi*d2vfdphi2/g;
    jac(3,2) = 0.75*pow(X[0],2)*dvfdtau/(X[3]*g)+X[3]*pow(X[4],2)*kPhi*dvfdtau/(6*g)+(3./16)*pow(X[4],2)*dkPhi*dvfdtau/g+(3./8)*pow(X[0],2)*d2vfdphidtau/g+(3./8)*pow(X[4],2)*kPhi*d2vfdphidtau/g;
    jac(3,3) = -5 -9/pow(X[3],2)+(4./3)*pow(X[3],2)-0.75*pow(X[0]/X[3],2)*vf/g+pow(X[4],2)*kPhi*vf/(6*g)+0.75*pow(X[0]/X[3],2)*vg;
    jac(3,4) = -0.75*X[4]*kPhi*vf/(X[3]*pow(g,3))-X[3]*pow(X[4],3)*pow(kPhi/X[0],2)*vf/(6*pow(g,3))+X[3]*X[4]*kPhi*vf/(3*g);
    jac(3,4) += -(3./16)*pow(X[4],3)*kPhi*vf*dkPhi/(pow(X[0],2)*pow(g,3))+(3./8)*X[4]*vf*dkPhi/g-(3./8)*X[4]*kPhi*dvfdphi/pow(g,3)-(3./8)*pow(X[4],3)*pow(kPhi/X[0],2)*dvfdphi/pow(g,3)+0.75*X[4]*kPhi*dvfdphi/g;
    jac(3,5) = 0; jac(3,6) = 0;
    jac(4,0) = -4*X[0]*X[2]/kPhi+8*pow(X[4]/X[0],3)*kPhi+pow(X[4],5)*pow(kPhi,2)*vf/(6*pow(X[0]*g,3))+X[3]*pow(X[4]/X[0],3)*dkPhi+2*X[3]*pow(X[4]/X[0],3)*kPhi*dvf0/vf0;
    jac(4,1) = -4*pow(X[4],3)*dkPhi/pow(X[0],2)+2*pow(X[0]/kPhi,2)*X[2]*dkPhi-pow(X[4],5)*kPhi*vf*dkPhi/(12*pow(X[0],2)*pow(g,3))+pow(X[4],3)*vf*dkPhi/(6*g)+X[3]*X[4]*pow(dkPhi/kPhi,2)-X[3]*pow(X[4],3)*dkPhi*dvf0/(pow(X[0],2)*vf0);
    jac(4,1) += X[3]*X[4]*pow(dvf0/vf0,2)+X[3]*pow(X[4],3)*kPhi*pow(dvf0/(X[0]*vf0),2)-0.5*X[3]*pow(X[4],3)*d2kPhi/pow(X[0],2)-X[3]*X[4]*d2kPhi/kPhi-X[3]*X[4]*d2vf0/vf0-X[3]*pow(X[4],3)*kPhi*d2vf0/(pow(X[0],2)*vf0)+pow(X[4],3)*kPhi*dvfdphi/(6*g);
    jac(4,2) = -2*pow(X[4],2)-2*pow(X[0],2)/kPhi + pow(X[4],3)*kPhi*dvfdtau/(6*g);
    jac(4,3) = (8./9)*X[3]*X[4]-pow(X[4],3)*dkPhi/(2*X[0]*X[0])-X[4]*dkPhi/kPhi - X[4]*dvf0/vf0-pow(X[4],3)*kPhi*dvf0/(pow(X[0],2)*vf0);
    jac(4,4) = -4+(4./9)*pow(X[3],2)-4*X[2]*X[4]-12*pow(X[4]/X[0],2)*kPhi-pow(X[4],4)*pow(kPhi/X[0],2)*vf/(6*pow(g,3));
    jac(4,4) += pow(X[4],2)*kPhi*vf/(2*g)-1.5*X[3]*pow(X[4]/X[0],2)*dkPhi-X[3]*dkPhi/kPhi-X[3]*dvf0/vf0-3*X[3]*pow(X[4]/X[0],2)*kPhi*dvf0/vf0;
    jac(4,5) = 0; jac(4,6) = 0;
    jac(5,0) = exp(-A); jac(5,1) = 0; jac(5,2) = 0; jac(5,3) = 0; jac(5,4) = 0; jac(5,5) = 0; jac(5,6) = 0;
    jac(6,0) = exp(-A)/g; jac(6,1) = exp(-A)*dkPhi*pow(X[4],2)/(2*X[0]*g); jac(6,2) = 0; jac(6,3) = 0; 
    jac(6,4) = exp(-A)*kPhi*X[4]/(X[0]*g); jac(6,5) = 0; jac(6,6) = 0;
    dfdt[0] = 0; dfdt[1] = 0; dfdt[2] = 0; dfdt[3] = 0; dfdt[4] = 0; dfdt[5] = -X[0]*exp(-A);
    dfdt[6] = - g*X[0]*exp(-A);
}

void HVQCD::eomCoupled(const state &X , state &dXdA , const double A)
{
    // X = (q, phi, tau, dphi, dtau, z, u)
    dXdA[0] = dqCoupled(X[0], X[1], X[2], X[3], X[4]);
    dXdA[1] = X[3];
    dXdA[2] = X[4];
    dXdA[3] = d2PhiCoupled(X[0], X[1], X[2], X[3], X[4]);
    dXdA[4] = d2tauCoupled(X[0], X[1], X[2], X[3], X[4]);
    dXdA[5] = dzdA(X[0], A);
    dXdA[6] = dudA(X[0], X[1], dXdA[2], A);
}

void HVQCD::observerCoupled(const state &X , const double A)
{
    // X = (q, Phi, tau, dPhi/dA, dtau/dA, z, u)
    double dqc = dqCoupled(X[0], X[1], X[2], X[3], X[4]);
    double d2phic = d2PhiCoupled(X[0], X[1], X[2], X[3], X[4]);
    double d2tc = d2tauCoupled(X[0], X[1], X[2], X[3], X[4]);
    qs.push_back(X[0]);
    Phis.push_back(X[1]);
    taus.push_back(X[2]);
    dqs.push_back(dqc);
    dPhis.push_back(X[3]);
    dtaus.push_back(X[4]);
    d2qs.push_back(d2qCoupled(X[0], X[1], X[2], dqc, X[3], X[4], d2tc));
    d2Phis.push_back(d2phic);
    d2taus.push_back(d2tc);
    d3taus.push_back(d3tauCoupled(X[0], X[1], X[2], dqc, X[3], X[4], d2phic, d2tc));
    As.push_back(A);
    zs.push_back(X[5]);
    us.push_back(X[6]);
}

double HVQCD::dqUV(const double q, const double lambda , const double dlambda)
{
    // Returns dq/dA from the UV EOMs of q, lambda and tau
    return (4.0/9) * q * pow(dlambda/lambda,2);
}

double HVQCD::d2qUV(const double q, const double lambda , const double dq, const double dlambda, const double d2lambda)
{
    // Returns d2q/dA2 from the UV EOMs of q, lambda and tau
    return (4.0/9)*dq*pow(dlambda/lambda,2)-(8./9)*q*pow(dlambda/lambda,3)+(8./9)*q*dlambda*d2lambda/pow(lambda,2);
}

double HVQCD::d2lambdaUV(const double q, const double lambda , const double dq, const double dlambda)
{
    // Returns d2lambda/dA2 from the UV EOMs of q, lambda and tau
   double vgl = Vgl(lambda);
   double dvgl = dVgldlambda(lambda);
   double vfl = Vfl(lambda, 0);
   double dvfdl = dVfldlambda(lambda, 0);
   double ans = -(3./8) * pow(q*lambda,2) * dvgl + 9 * pow(lambda,2)/dlambda + (3./4) * pow(q*lambda,2) * (vfl-vgl)/dlambda ;
   ans = ans - 5 * dlambda + dq * dlambda / q + pow(dlambda,2)/lambda + (3./8) * pow(q*lambda,2) * dvfdl;
   return ans;
}

double HVQCD::d2taunUV(const double q, const double lambda, const double tau, const double dq, const double dlambda, const double dtau)
{
    // Returns d2taun/dA2 in the UV regime
    double vf0l = Vf0l(lambda);
    double dvf0dl = dVf0dlambda(lambda);
    double kl = klambda(lambda);
    double dkdl = dkdlambda(lambda);
    double ans = 3 * tau - 2 * q * q * tau / kl - tau * dq / q + tau * dvf0dl * dlambda / vf0l;
    ans = ans - 2 * dtau  + dq * dtau / q - dvf0dl * dlambda * dtau / vf0l + tau * dlambda * dkdl / kl;
    ans = ans - dlambda * dtau * dkdl / kl;
    return ans; 
}

double HVQCD::d3taunUV(const double q, const double lambda, const double tau, const double dq, const double dlambda, const double dtau,
                       const double d2q, const double d2lambda, const double d2tau)
{
    // Returns d3taun/dA3 in the UV regime
    double kl = klambda(lambda);
    double dkl = dkdlambda(lambda);
    double d2kl = d2kdlambda2(lambda);
    double vf0 = Vf0l(lambda);
    double dvf0dl = dVf0dlambda(lambda);
    double d2vf0dl2 = d2Vf0dlambda2(lambda);
    double ans = -4*q*tau*dq/kl+tau*pow(dq/q,2)+2*pow(q/kl,2)*tau*dkl*dlambda-tau*pow(dvf0dl*dlambda/vf0,2);
    ans += -tau*pow(dkl/kl*dlambda,2)+3*dtau-2*q*q*dtau/kl-dq*dtau/q-pow(dq/q,2)*dtau;
    ans += dvf0dl*dlambda*dtau/vf0+dkl*dlambda*dtau/kl+pow(dvf0dl*dlambda/vf0,2)*dtau;
    ans += pow(dkl*dlambda/kl,2)*dtau-tau*d2q/q+dtau*d2q/q+tau*pow(dlambda,2)*d2vf0dl2/vf0;
    ans += -pow(dlambda,2)*dtau*d2vf0dl2/vf0+tau*pow(dlambda,2)*d2kl/kl-pow(dlambda,2)*dtau*d2kl/kl;
    ans += tau*dvf0dl*d2lambda/vf0 + tau*dkl*d2lambda/kl - dvf0dl*dtau*d2lambda/vf0;
    ans += -dkl*dtau*d2lambda/kl -2*d2tau+dq*d2tau/q-dvf0dl*dlambda*d2tau/vf0-dkl*dlambda*d2tau/kl;
    return ans;
}

void HVQCD::eomUV(const state &X, state &dXdA, const double A)
{
    // X = (q, lambda, tau, dlambda, dtau)
    dXdA[0] = dqUV(X[0], X[1] , X[3]);
    dXdA[1] = X[3];
    dXdA[2] = X[4];
    dXdA[3] = d2lambdaUV(X[0], X[1], dXdA[0], X[3]);
    dXdA[4] = d2taunUV(X[0], X[1], X[2], dXdA[0], X[3], X[4]);
} 

void HVQCD::observerUV(const state &X , double A)
{
    // X = q, lambda, taunUV, dlambda, dtaunUV
    qs.push_back(X[0]);
    Phis.push_back(log(X[1]));                  // Phi = log(lambda)
    lUVs.push_back(X[1]);
    taus.push_back(X[2]*exp(-A));               // tau = exp(-A) taun
    tauns.push_back(X[2]);
    double dq = dqUV(X[0], X[1] , X[3]);
    dqs.push_back(dq);
    dPhis.push_back(X[3]/X[1]);                 // dPhi = dlambda / lambda
    dlUVs.push_back(X[3]);
    dtaus.push_back(exp(-A)*(X[4]-X[2]));       // dtau = exp(-A)(dtaun -taun)
    dtauns.push_back(X[4]);
    double d2lambda = d2lambdaUV(X[0], X[1] , dq, X[3]);
    double d2q = d2qUV(X[0], X[1], dq, X[3], d2lambda);
    d2qs.push_back(d2q);
    d2Phis.push_back(d2lambda/X[1] - pow(X[3]/X[1],2));
    double d2taun = d2taunUV(X[0], X[1], X[2], dq, X[3], X[4]);
    d2taus.push_back(exp(-A)*(X[2]-2*X[4]+d2taun));    // d2tau = exp(-A)(taun-2dtaun+d2taun)
    double d3taun = d3taunUV(X[0], X[1], X[2], dq, X[3], X[4], d2q, d2lambda, d2taun);
    d3taus.push_back(exp(-A)*(-X[2]+3*X[4]-3*d2taun+d3taun));
    As.push_back(A);
    AUVs.push_back(A);
}

std::vector<double> HVQCD::computeVectorMesonPotential()
{
    // Returns the potential of the Flavour Non-Singlet Vector Mesons
    std::vector<double> V(As.size());
    for(int i = 0; i < As.size(); i++)
    {
        double e2A = exp(2*As[i]);
        double g = G(qs[i], Phis[i], dtaus[i]);
        double dg = dG(qs[i], Phis[i], dqs[i], dPhis[i], dtaus[i], d2taus[i]);
        double wPhi = w(Phis[i]);
        double dwPhi = dwdPhi(Phis[i]);
        double d2wPhi = d2wdPhi2(Phis[i]);
        double vf0 = Vf0(Phis[i]);
        double dvf0 = dVf0dPhi(Phis[i]);
        double d2vf0 = d2Vf0dPhi2(Phis[i]);
        V[i] = 0.75 - 0.5 * dg/g - 0.5*dqs[i]/qs[i] -2*taus[i]*dtaus[i]+taus[i]*dg*dtaus[i]/g+taus[i]*dqs[i]*dtaus[i]/qs[i];
        V[i] = V[i] - pow(dtaus[i],2) + pow(taus[i]*dtaus[i],2) + dvf0*dPhis[i]/vf0 - 0.5*dg*dvf0*dPhis[i]/(g*vf0)-0.5*dqs[i]*dvf0*dPhis[i]/(qs[i]*vf0);
        V[i] += 2*dwPhi*dPhis[i]/wPhi - dg*dwPhi*dPhis[i]/(g*wPhi) - dqs[i]*dwPhi*dPhis[i]/(qs[i]*wPhi)-taus[i]*dvf0*dtaus[i]*dPhis[i]/vf0;
        V[i] += -2*taus[i]*dwPhi*dtaus[i]*dPhis[i]/wPhi - 0.25*pow(dvf0*dPhis[i]/vf0,2)+dvf0*dwPhi*pow(dPhis[i],2)/(vf0*wPhi);
        V[i] += 0.5*pow(dPhis[i],2)*d2vf0/vf0 + pow(dPhis[i],2)*d2wPhi/wPhi-taus[i]*d2taus[i]+0.5*dvf0*d2Phis[i]/vf0 + dwPhi*d2Phis[i]/wPhi;
        V[i] = e2A*V[i]/pow(g*qs[i],2);
    }
    return V;
}

std::vector<double> HVQCD::computeAxialVectorMesonNonSingletPotential(const std::vector<double> &VVectorMeson)
{
    // Computation of the flavour Non-singlet axial mesons' potential
    std::vector<double> VAxialMesonsNonSinglet(As.size());
    for(int i = 0; i < As.size(); i++)
    {
        double wPhi = w(Phis[i]);
        double kPhi = k(Phis[i]);
        VAxialMesonsNonSinglet[i] = VVectorMeson[i] + pow(2*taus[i]*exp(As[i])/wPhi,2)*kPhi ;
    }
    return VAxialMesonsNonSinglet;
}

std::vector<double> HVQCD::computePseudoScalarMesonPotential()
{
    // Computes the Pseudoscalar Meson potential
    std::vector<double> V(As.size());
    for(int i = 0; i < As.size(); i++)
    {
        double e2A = exp(2*As[i]);
        double g = G(qs[i], Phis[i], dtaus[i]);
        double dg = dG(qs[i], Phis[i], dqs[i], dPhis[i], dtaus[i], d2taus[i]);
        double kPhi = k(Phis[i]);
        double dkPhi = dkdPhi(Phis[i]);
        double d2kPhi = d2kdPhi2(Phis[i]);
        double wPhi = w(Phis[i]);
        double vf0 = Vf0(Phis[i]);
        double dvf0 = dVf0dPhi(Phis[i]);
        double d2vf0 = d2Vf0dPhi2(Phis[i]);
        V[i] = 0.75+1.5*dg/g+1.5*dqs[i]/qs[i]+2*dtaus[i]/taus[i]-2*taus[i]*dtaus[i]+dg*dtaus[i]/(g*taus[i])-taus[i]*dg*dtaus[i]/g;
        V[i] += dqs[i]*dtaus[i]/(qs[i]*taus[i])-taus[i]*dqs[i]*dtaus[i]/qs[i]-pow(dtaus[i],2)+2*pow(dtaus[i]/taus[i],2)+pow(taus[i]*dtaus[i],2)+dkPhi*dPhis[i]/kPhi;
        V[i] += dg*dkPhi*dPhis[i]/(2*g*kPhi)+dkPhi*dqs[i]*dPhis[i]/(2*kPhi*qs[i])+dvf0*dPhis[i]/vf0+dg*dvf0*dPhis[i]/(2*g*vf0);
        V[i] += dqs[i]*dvf0*dPhis[i]/(2*qs[i]*vf0)+dkPhi*dtaus[i]*dPhis[i]/(kPhi*taus[i])-taus[i]*dkPhi*dtaus[i]*dPhis[i]/kPhi;
        V[i] += dvf0*dtaus[i]*dPhis[i]/(vf0*taus[i])-taus[i]*dvf0*dtaus[i]*dPhis[i]/vf0+0.75*pow(dkPhi*dPhis[i]/kPhi,2);
        V[i] += dkPhi*dvf0*pow(dPhis[i],2)/(2*kPhi*vf0)+0.75*pow(dvf0*dPhis[i]/vf0,2)-pow(dPhis[i],2)*d2kPhi/(2*kPhi)-pow(dPhis[i],2)*d2vf0/(2*vf0);
        V[i] += -d2taus[i]/taus[i] + taus[i]*d2taus[i]-dkPhi*d2Phis[i]/(2*kPhi)-dvf0*d2Phis[i]/(2*vf0);
        V[i] = e2A*V[i]/pow(g*qs[i],2) + pow(2*taus[i]*exp(As[i])/wPhi,2)*kPhi;
    }
    return V;   
}

std::vector<double> HVQCD::computeScalarMesonPotential()
{
    std::vector<double> V(As.size());
    for(int i = 0; i < As.size(); i++)
    {
        double e2A = exp(2*As[i]);
        double g = G(qs[i], Phis[i], dtaus[i]);
        double dg = dG(qs[i], Phis[i], dqs[i], dPhis[i], dtaus[i], d2taus[i]);
        double d2g = d2G(qs[i], Phis[i], dqs[i], dPhis[i], dtaus[i], d2qs[i], d2Phis[i], d2taus[i], d3taus[i]);
        double kPhi = k(Phis[i]);
        double dkPhi = dkdPhi(Phis[i]);
        double d2kPhi = d2kdPhi2(Phis[i]);
        double wPhi = w(Phis[i]);
        double vf0 = Vf0(Phis[i]);
        double dvf0 = dVf0dPhi(Phis[i]);
        double d2vf0 = d2Vf0dPhi2(Phis[i]);
        V[i] = 3.75 - 5.5*dg/g + 3*pow(dg/g,2) -1.5*dqs[i]/qs[i] + dg*dqs[i]/(g*qs[i])-4*taus[i]*dtaus[i];
        V[i] += 3*taus[i]*dg*dtaus[i]/g + taus[i]*dqs[i]*dtaus[i]/qs[i]-pow(dtaus[i],2)+pow(taus[i]*dtaus[i],2)+2*dkPhi*dPhis[i]/kPhi;
        V[i] += -1.5*dg*dkPhi*dPhis[i]/(g*kPhi)-dkPhi*dqs[i]*dPhis[i]/(2*kPhi*qs[i])+2*dvf0*dPhis[i]/vf0-1.5*dg*dvf0*dPhis[i]/(g*vf0);
        V[i] += -dqs[i]*dvf0*dPhis[i]/(2*qs[i]*vf0)-taus[i]*dkPhi*dtaus[i]*dPhis[i]/kPhi-taus[i]*dvf0*dtaus[i]*dPhis[i]/vf0;
        V[i] += -0.25*pow(dkPhi*dPhis[i]/kPhi,2) + dkPhi*dvf0*pow(dPhis[i],2)/(2*kPhi*vf0)-0.25*pow(dvf0*dPhis[i]/vf0,2)-d2g/g;
        V[i] += pow(dPhis[i],2)*d2kPhi/(2*kPhi)+pow(dPhis[i],2)*d2vf0/(2*vf0)-taus[i]*d2taus[i]+dkPhi*d2Phis[i]/(2*kPhi)+dvf0*d2Phis[i]/(2*vf0);
        V[i] = e2A*V[i]/pow(qs[i]*g,2) -2*e2A/kPhi;
    }
    return V;
}

std::vector<double> HVQCD::computeAxialVectorMesonSingletPotential(const std::vector<double> &VAxialVectorMeson)
{
    std::vector<double> V(As.size());
    for(int i = 0; i < As.size(); i++)
    {
        double e2A = exp(2*As[i]);
        double l = std::exp(Phis[i]);
        double z = Z(l);
        double vf0 = Vf0(Phis[i]);
        double g = G(qs[i], Phis[i], dtaus[i]);
        double wPhi = w(Phis[i]);
        V[i] = VAxialVectorMeson[i] + 4.0 * xf * e2A * z * std::exp(-taus[i]*taus[i])/(vf0 * g * wPhi * wPhi);
    }
    return V;
}

void HVQCD::finalizeBackground()
{
    // Selects the relevant values to compute the potentials later
    std::vector<double> A, z, u, q, Phi, tau, dq, dPhi, dtau, d2q, d2Phi, d2tau, d3tau;
    for(int i = 0; i < As.size(); i++)
    {
        // For the spectrum we are only interested in -50 < A < 20
        if ((zs[i]-zs.back()> 1e-6) && (Phis[i] < 120))
        {
            A.push_back(As[i]); z.push_back(zs[i]-zs.back()); u.push_back(us[i]-us.back());
            q.push_back(qs[i]); Phi.push_back(Phis[i]); tau.push_back(taus[i]);
            dq.push_back(dqs[i]); dPhi.push_back(dPhis[i]); dtau.push_back(dtaus[i]);
            d2q.push_back(d2qs[i]); d2Phi.push_back(d2Phis[i]); d2tau.push_back(d2taus[i]);
            d3tau.push_back(d3taus[i]);
        }
    }
    As = A; zs = z; us = u; qs = q; Phis = Phi, taus = tau; dqs = dq; dPhis = dPhi; dtaus = dtau;
    d2qs = d2q; d2Phis = d2Phi; d2taus = d2tau;
    d3taus = d3tau;
}

void HVQCD::solve()
{
   // Solves the Holographic QCD model
   // Clear the vector containers
   qs.clear(); Phis.clear(); taus.clear(); dqs.clear(); dPhis.clear(); dtaus.clear();
   d2qs.clear(); d2Phis.clear(); d2taus.clear(); d3taus.clear(); As.clear(); zs.clear(); us.clear();
   AYM2.clear(); zYM2.clear(); uYM2.clear(); qYM2.clear(); PhiYM2.clear(); tauYM2.clear();
   dqYM2.clear(); dPhiYM2.clear(); dtauYM2.clear(); d2qYM2.clear(); d2PhiYM2.clear();
   d2tauYM2.clear(); d3tauYM2.clear();
   lUVs.clear(); tauns.clear(); dlUVs.clear(); dtauns.clear(); AUVs.clear();
   // Setting up IR boundary conditions
   double Air = -150, AUVYM = 50.0, h = 0.01, AUVc = 100, AUVf = 1000;
   std::function<double(double)> func = [this, &Air] (double z) { return this->AIR(z) - Air ;} ;
   double zIRYM = zbrent(func, 0.1, 100.0, 1e-9, true), uIR = zIRYM;
   double tcoeff = (12 - xf * W0) * kIR / (6.0 * VgIR) ;
   double qir = exp(Air) / dAIR(zIRYM), Phiir = PhiIR(zIRYM), tauIR = tau0 * pow(zIRYM, tcoeff);
   double dqir = dqYM(qir, Phiir), dPhiir = dPhiYM(qir, Phiir), dtauIR = tcoeff * tau0 * pow(zIRYM, tcoeff - 1.0) * exp(-Air) * qir;
   double d2Phiir = d2PhiYM(qir, Phiir), d2tauIR = exp(-2*Air)*tcoeff*tau0*pow(zIRYM,tcoeff-2)*((tcoeff-1)*qir*qir-exp(Air)*zIRYM*(qir-dqir));
   double d2qir = d2qYM(qir, Phiir);
   double d3tauIR = d3tauYM(qir, Phiir, tauIR, dqir, dPhiir, dtauIR, d2qir, d2Phiir,d2tauIR);
   std::function<double(double)> tcorr = [this] (double l)
   {
        return (-88+16*this->xf+27*this->sc*this->kU1)*log(24*pow(M_PI,2)/((11-2*this->xf)*l))/(-66+12*this->xf);
   };
   // Get Yang Mills profile of q and Phi
   double tcut = 1000, Vfcut = 1e-8;
   state X(4);
   X <<= qir, Phiir, tauIR, zIRYM;
   typedef boost::numeric::odeint::result_of::make_dense_output<boost::numeric::odeint::runge_kutta_dopri5<state> >::type dense_stepper;
   if (X[2] > tcut)
   {
        double A = Air, q = qir, phi = Phiir, tau = tauIR, z = zIRYM, u = uIR;
        double dq = dqir, dphi = dPhiir, dtau =  dtauYangMills1(qir, Phiir, tauIR, dPhiir);
        double d2q = d2qir, d2phi = d2Phiir, d2tau = d2tauIR, d3tau = d3tauIR;
        dense_stepper stepper_IR_CUT = make_dense_output(1.0e-12, 1.0e-12, boost::numeric::odeint::runge_kutta_dopri5< state >());
        stepper_IR_CUT.initialize(X, A, h);
        auto eomfunYM1 = [this] (const state &Y, state &dYdA, double A) {this->eomYangMills1(Y, dYdA, A);};
        while (tau > tcut)
        {
            // Add the fields to the containers
            As.push_back(A); us.push_back(u); zs.push_back(z); qs.push_back(q); Phis.push_back(phi);
            taus.push_back(tau); dqs.push_back(dq); dPhis.push_back(dphi); dtaus.push_back(dtau);
            d2qs.push_back(d2q); d2Phis.push_back(d2phi); d2taus.push_back(d2tau); d3taus.push_back(d3tau);
            stepper_IR_CUT.do_step(eomfunYM1);
            A = stepper_IR_CUT.current_time();
            X = stepper_IR_CUT.current_state();
            q = X[0]; phi = X[1]; tau = X[2]; z = X[3]; u = X[4];
            dq = dqYM(q, phi); dphi = dPhiYM(q, phi);
            dtau = dtauYangMills1(q, phi, X(2), dphi);
            d2q = d2qYM(q, phi); d2phi = d2PhiYM(q, phi);
            d2tau = d2tauYangMills1(q, phi, tau, dq, dphi, d2phi);
            d3tau = d3tauYM(q, phi, tau, dq, dphi, dtau, d2q, d2phi, d2tau);
        }
   }
   // Now we want to solve the tachyon from A until
   // Vf0(PhiYM)exp(-tau^2) == Vfcut Vg(PhiYM)
   // We first compute qYM, PhiYM and tauYM2 simultaneously
   // Then we solve the above equation
   // First we define the boundary conditions
   double AUV1; state XtauYM2(6);
   if(As.size() != 0)
   {
       AUV1 = As.back();
       XtauYM2 <<= qs.back(), Phis.back(), taus.back(), dtaus.back(), zs.back(), us.back();
   }
   else
   {
       AUV1 = Air;
       XtauYM2 <<= qir, Phiir, tauIR, dtauIR, zIRYM, uIR;
   }
   auto eomfunYM2 = [this] (const state &Y, state &dYdA, double A) {this->eomYangMills2(Y, dYdA, A);};
   auto obsfunYM2 = [this] (const state &Y, double A) {this->observerYangMills2(Y, A);};
   // Solve the EOMs
   dense_stepper stepperYM2 = make_dense_output(1.0e-12, 1.0e-12, boost::numeric::odeint::runge_kutta_dopri5< state >());
   integrate_const(stepperYM2, eomfunYM2, XtauYM2, AUV1, AUVYM, h, obsfunYM2);
   // Get the profiles
   Spline_Interp<double> qYM2Profile(AYM2, qYM2, dqYM2.front(), dqYM2.back());
   Spline_Interp<double> PhiYM2Profile(AYM2, PhiYM2, dPhiYM2.front(), dPhiYM2.back());
   Spline_Interp<double> tauYM2Profile(AYM2, tauYM2, dtauYM2.front(), dtauYM2.back());
   Spline_Interp<double> zYM2Profile(AYM2, zYM2, qYM2.front()*exp(-AYM2.front()),qYM2.back()*exp(-AYM2.back()));
   double duIR = dudA(qYM2.front(), PhiYM2.front(), dtauYM2.front(), AYM2.front());
   double duUV = dudA(qYM2.back(), PhiYM2.back(), dtauYM2.back(), AYM2.back());
   Spline_Interp<double> uYM2Profile(AYM2, uYM2, duIR, duUV);
   // We now compute AUV2 such that Vf0(PhiYM)exp(-tau^2) == Vfcut Vg(PhiYM)
   std::function<double(double)> func2 = [this, &Vfcut, &PhiYM2Profile, &tauYM2Profile] (double A) 
   { 
    double phi = PhiYM2Profile.interp(A);
    double tau = tauYM2Profile.interp(A);
    double ans = this->Vf0(phi) * exp(-tau * tau)  - Vfcut * this->Vg(phi);
    return ans;
   };
   double AUV2 = zbrent(func2, AUV1, AUVYM, 1e-9, true);
   double PhiUV2 = PhiYM2Profile.interp(AUV2), tauUV2 = tauYM2Profile.interp(AUV2);
   bool potCond = fabs(Vf0(PhiUV2) * exp(-pow(tauUV2,2))  - Vfcut * Vg(PhiUV2)) > (0.01 * Vfcut * Vg(PhiUV2));
   if (AUV2 < AUV1 || AUV2 > AUVYM || potCond)
   {
       As.insert(As.end(), AYM2.begin(), AYM2.end());
       zs.insert(zs.end(), zYM2.begin(), zYM2.end());
       us.insert(us.end(), uYM2.begin(), uYM2.end());
       qs.insert(qs.end(), qYM2.begin(), qYM2.end());
       dqs.insert(dqs.end(),dqYM2.begin(), dqYM2.end());
       Phis.insert(Phis.end(), PhiYM2.begin(), PhiYM2.end());
       dPhis.insert(dPhis.end(), dPhiYM2.begin(), dPhiYM2.end());
       taus.insert(taus.end(), tauYM2.begin(), tauYM2.end());
       dtaus.insert(dtaus.end(), dtauYM2.begin(), dtauYM2.end());
       d2qs.insert(d2qs.end(), d2qYM2.begin(), d2qYM2.end());
       d2Phis.insert(d2Phis.end(), d2PhiYM2.begin(), d2PhiYM2.end());
       d2taus.insert(d2taus.end(), d2tauYM2.begin(), d2tauYM2.end());
       d3taus.insert(d3taus.end(), d3tauYM2.begin(), d3tauYM2.end());
       finalizeBackground();
       return;
   }
   // Append the values and find iYM2 that will allow us
   // to start to integrate the coupled EOMs
   int iYM2;
   for(int i = 0; i < AYM2.size(); i++)
   {
       if (AYM2[i] < AUV2)
       {
           As.push_back(AYM2[i]); zs.push_back(zYM2[i]); us.push_back(uYM2[i]); qs.push_back(qYM2[i]);
           Phis.push_back(PhiYM2[i]); taus.push_back(tauYM2[i]); dqs.push_back(dqYM2[i]); dPhis.push_back(dPhiYM2[i]);
           dtaus.push_back(dtauYM2[i]); d2qs.push_back(d2qYM2[i]); d2Phis.push_back(d2PhiYM2[i]); d2taus.push_back(d2tauYM2[i]);
           d3taus.push_back(d3tauYM2[i]);
           iYM2 = i;
       }
   }
   // We now solve the coupled EOMs
   /* First we set up the BCs. For that we need to determine
      the q, Phi, tau, dPhi and dtau at AUV2
      For that we use profile computed previously
   */
   // Boundary conditions
   state Xcoupled(7);
   Xcoupled <<= qYM2Profile.interp(AUV2), PhiYM2Profile.interp(AUV2), tauYM2Profile.interp(AUV2), PhiYM2Profile.der1(AUV2), tauYM2Profile.der1(AUV2), zYM2Profile.interp(AUV2), uYM2Profile.interp(AUV2);
   auto jacfun = [this] (const state &Y , matrix_type &jac , const double A, state &dfdt) {this->jacobian(Y, jac, A, dfdt);};
   auto eomcoupled = [this] (const state &Y, state &dYdA, const double A) {this->eomCoupled(Y, dYdA, A);};
   auto obscoupled = [this] (const state &Y, const double A) {this->observerCoupled(Y, A);};
   typedef boost::numeric::odeint::result_of::make_dense_output<boost::numeric::odeint::rosenbrock4<double> >::type dense_stepper_stiff;
   dense_stepper_stiff stepper_AUV2_AUVC = make_dense_output(1.0e-9, 1.0e-9, boost::numeric::odeint::rosenbrock4< double >());
   integrate_adaptive(stepper_AUV2_AUVC, std::make_pair(eomcoupled, jacfun), Xcoupled, AUV2, AUVc, h, obscoupled);
   double AUV3 = As.back();
   if( AUV3 < AUVc - 0.1 )
   {
       finalizeBackground();
       return ;
   }
   // Solving the UV EOMs
   // Boundary conditions
   state XUV(5);
   XUV <<= qs.back(), exp(Phis.back()), exp(As.back()) * taus.back(), dPhis.back() * exp(Phis.back()), exp(As.back()) * (taus.back() + dtaus.back()) ;
   auto uveom = [this] (const state &Y, state &dYdA, const double A) {this->eomUV(Y, dYdA, A);};
   auto obsUV = [this] (const state &Y, const double A) {this->observerUV(Y, A);};
   dense_stepper stepper_AUVC_AUVF = make_dense_output(1.0e-12, 1.0e-12, boost::numeric::odeint::runge_kutta_dopri5< state >());
   integrate_adaptive(stepper_AUVC_AUVF, uveom, XUV, AUV3, AUVf, h, obsUV);
   if( As.back() < AUVf - 0.1)
   {
       finalizeBackground();
       return ;
   }
   // Now we compute mq
   Spline_Interp<double> taunUV(AUVs, tauns, dtauns.front(), dtauns.back());
   Spline_Interp<double> lambdaUV(AUVs, lUVs, dlUVs.front(), dlUVs.back());
   double lUV = 1;
   mq = lambdaUV.interp(AUVf-10) * taunUV.interp(AUVf)*exp(-log(lUV) - tcorr(lambdaUV.interp(AUVf))) - lambdaUV.interp(AUVf) * taunUV.interp(AUVf - 10)*exp(-log(lUV) - tcorr(lambdaUV.interp(AUVf - 10)));
   mq = mq / (lambdaUV.interp(AUVf-10) - lambdaUV.interp(AUVf)) ;
   finalizeBackground();
   return ;
}

void HVQCD::computeSpectrum(std::string method)
{
    /*
        Computes the lowest 2++ glueball, first 6 Nonsinglet Vector Meson, first 5 Nonsinglet Axial Vector Meson,
        first 5 Nonsinglet Pseudoscalar Meson and first 2 Nonsinglet scalar meson masses of HVQCD.
        If HVQCD is not solved this function solves it first and then computes the spectrum. The results are then printed.
    */
    // Check if the theory is solved
    if (As.size() == 0) this->solve();
    // Compute the Schrodinger potential of the fluctuations
    std::vector<double> V2G, VVM, VAVM, VPSM, VSM, VSingletAVM;
    if (add_tensor_glueball) V2G = computeV2GPotential();
    VVM = computeVectorMesonPotential();
    VAVM = computeAxialVectorMesonNonSingletPotential(VVM);
    if (add_scalars)
    {
        VPSM = computePseudoScalarMesonPotential();
        VSM = computeScalarMesonPotential();
    }
    if (add_singlet_axial) VSingletAVM = computeAxialVectorMesonSingletPotential(VAVM);
    // Compute the masses
    std::vector<double> TGMasses, VMMasses, AVMMasses, PSMMasses, SMMasses, SingletAVMMasses;
    if (add_tensor_glueball) TGMasses = computeMasses(zs, V2G, 1, method);
    VMMasses = computeMasses(us, VVM, 6, method);
    AVMMasses = computeMasses(us, VAVM, 5, method);
    if(add_scalars)
    {
        PSMMasses = computeMasses(us, VPSM, 5, method);
        SMMasses = computeMasses(us, VSM, 2, method);
    }
    if(add_singlet_axial) SingletAVMMasses = computeMasses(us, VSingletAVM, 2, method);
    // Display the mass values
    // Tensor glueball ratio
    if(add_tensor_glueball)
    {
        std::cout << "TENSOR GLUEBALL SECTOR" << std::endl;
        std::cout << TGMasses[0] << std::endl;
        std::cout << std::endl;
    }
    std::cout << "VECTOR MESON NONSINGLET SECTOR" << std::endl;
    for(int i = 0; i < VMMasses.size(); i++) std::cout << VMMasses[i] << '\t';
    std::cout << std::endl << std::endl;
    std::cout << "AXIAL VECTOR MESON NONSINGLET SECTOR" << std::endl;
    for(int i = 0; i < AVMMasses.size(); i++) std::cout << AVMMasses[i] << '\t';
    std::cout << std::endl << std::endl;
    if(add_scalars)
    {
        std::cout << "PSEUDOSCALAR MESON NONSINGLET SECTOR" << std::endl;
        for(int i = 0; i < PSMMasses.size(); i++) std::cout << PSMMasses[i] << '\t' ;
        std::cout << std::endl << std::endl;
        std::cout << "SCALAR MESON NONSINGLET SECTOR" << std::endl;
        for(int i = 0; i < SMMasses.size(); i++) std::cout << SMMasses[i] << '\t';
        std::cout << std::endl << std::endl;
    }
    // In this model the Singlet and Nonsinglet Vector Mesons have the same mass
    if(add_singlet_vector)
    {
        std::cout << "VECTOR MESON SINGLET SECTOR" << std::endl;
        for(int i = 0; i < VMMasses.size(); i++) std::cout << VMMasses[i] << '\t';
        std::cout << std::endl << std::endl;
    }
    if(add_singlet_axial)
    {
        std::cout << "AXIAL VECTOR MESON SINGLET SECTOR" << std::endl;
        for(int i = 0; i < SingletAVMMasses.size(); i++) std::cout << SingletAVMMasses[i] << '\t';
        std::cout << std::endl << std::endl;
    }
}

double HVQCD::J()
{
    std::cout << "sc: " << sc << " ksc: " << ksc << " wsc: " << wsc << " W0: " << W0 << " w0: " << w0 << " kU1: " << kU1;
    std::cout << " wU1: " << wU1 <<  " VgIR: " << VgIR << " WIR: " << WIR << " kIR: " << kIR << " wIR: " << wIR << " W1: " << W1;
    std::cout << " k1: " << k1 << " w1: " << w1 << " xf: " << xf << " tau0: " << tau0 << " Za: " << Za << " ca: " << ca << std::endl;
    // Solve the background
    try
    {
        this->solve();
    }
    catch(...)
    {
        std::cout << "Bad background:" << std::endl;
        std::cout << "erms: " << 1e99 << std::endl;
        return 1e99;
    }
    // Compute the potentials
    std::vector<double> V2G, VVM, VAVM, VPSM, VSM, VSingletAVM;
    try
    {
        if(add_tensor_glueball) V2G = computeV2GPotential();
        VVM = computeVectorMesonPotential();
        VAVM = computeAxialVectorMesonNonSingletPotential(VVM);
        if (add_scalars)
        {
            VPSM = computePseudoScalarMesonPotential();
            VSM = computeScalarMesonPotential();
        }
        if (add_singlet_axial) VSingletAVM = computeAxialVectorMesonSingletPotential(VAVM);
    }
    catch(...)
    {
        std::cout << "Potentials not computed, check these values" << std::endl;
        std::cout << "erms: " << 1e99 << std::endl;
        return 1e99;
    }
    // Compute the masses
    std::vector<double> TGMasses, VMMasses, AVMMasses, PSMMasses, SMMasses, SingletAVMMasses;
    try
    {
        if (add_tensor_glueball) TGMasses = computeMasses(zs, V2G, 1, "cheb");
        VMMasses = computeMasses(us, VVM, 6, "cheb");
        AVMMasses = computeMasses(us, VAVM, 5, "cheb");
        if (add_scalars)
        {
            PSMMasses = computeMasses(us, VPSM, 5, "cheb");
            SMMasses = computeMasses(us,VSM, 2, "cheb");
        }
        if(add_singlet_axial)
        {
            SingletAVMMasses = computeMasses(us, VSingletAVM, 2, "cheb");
            if(SingletAVMMasses.size() == 0)
            {
                std::cout << "Only negative values for Singlet AVM masses" << std::endl;
                throw "error";
            }
        }
    }
    catch(...)
    {
        std::cout << "Unable to compute spectrum. Check these values" << std::endl;
        std::cout << "erms: " << 1e99 << std::endl;
        return 1e99;
    }
    double erms = 0;
    // Chi2 contribution due to the tensor glueball ratio
    if (add_tensor_glueball) erms += std::pow((TGMasses[0]/VMMasses[0] - RTG_rho[0])/RTG_rho[0], 2);
    // Contributions from the Meson ratios
    // Contributions from the ratios m_{\rho_n}/m_{\rho}
    for(int i = 0; i < Rrho_rho.size(); i++) erms += std::pow((VMMasses[i+1]/VMMasses[0]-Rrho_rho[i])/Rrho_rho[i],2);
    // Contributions from the ratios m_{\a1_n}/m_{\rho}
    for(int i = 0; i < Ra1_rho.size(); i++) erms += std::pow((AVMMasses[i]/VMMasses[0]-Ra1_rho[i])/Ra1_rho[i],2);
    if(add_scalars)
    {
        // Contributions from the ratios m_{\pi_n} / m_{\rho}
        for(int i = 0; i < Rpi_rho.size(); i++) erms += std::pow((PSMMasses[i]/VMMasses[0]-Rpi_rho[i])/Rpi_rho[i],2);
        // Contributions from the ratios m_{a0_n} / m_{\rho}
        for(int i = 0; i < Ra0_rho.size(); i++) erms += std::pow((SMMasses[i]/VMMasses[0]- Ra0_rho[i])/Ra0_rho[i], 2);
    }
    // Singlet vector and axial vector meson sector
    if (add_singlet_vector) for(int i = 0; i < Romega_rho.size(); i++) erms += std::pow((VMMasses[i]/VMMasses[0]-Romega_rho[i])/Romega_rho[i],2);
    if (add_singlet_axial) for(int i = 0; i < Rf1_rho.size(); i++) erms += std::pow((SingletAVMMasses[i]/VMMasses[0]-Rf1_rho[i])/Rf1_rho[i],2);
    // We want to penalize backgrounds with mq != 0
    if (add_mass) erms += std::pow(mq,2.0) ;
    int nRatios = 0;
    if(add_tensor_glueball) nRatios += 1;
    if(add_scalars) nRatios += 17;
    else nRatios += 10;
    if(add_singlet_vector) nRatios += 3;
    if(add_singlet_axial) nRatios += 2;

    erms = sqrt(erms)/nRatios;

    if (std::isnan(erms))
    {
        erms = 1e99;
        std::cout << "None erms" << std::endl;
        std::cout << "erms: " << erms << std::endl;
    }
    else
    {
        std::cout << "mq: " << mq << " erms: " << erms << std::endl;
    }
    
    return erms;
}

void HVQCD::printQuarkMass()
{
    // Prints the quark mass mq
    if (As.size()!= 0) std::cout << "mq: " << mq << std::endl;
    else std::cout << "Solve the background first" << std::endl;
}

void HVQCD::saveBackgroundFields(std::string path)
{
    // Saves the background fields in a file at path
    std::ofstream myfile;
    myfile.open(path);
    /*myfile << "sc" << '\t' << "ksc" << '\t' << "wsc" << '\t';
    myfile << "W0" << '\t' << "w0" << '\t' << "kU1" << '\t' << "wU1" << '\t';
    myfile << "VgIR" << '\t' << "WIR" << '\t' << "kIR" << '\t' << "wIR" << '\t';
    myfile << "W1" << '\t' << "k1" << '\t' << "w1" << '\t' << "xf" << '\t' << "tau0" << std::endl;
    myfile << sc << '\t' << ksc << '\t' << wsc << '\t';
    myfile << W0 << '\t' << w0 << '\t' << kU1 << '\t' << wU1 << '\t';
    myfile << VgIR << '\t' << WIR << '\t' << kIR << '\t' << wIR << '\t';
    myfile << W1 << '\t' << k1 << '\t' << w1 << '\t' << xf << '\t' << tau0 << std::endl;
    myfile << std::endl;*/
    myfile << "A" << '\t' << "q" << '\t' << "Phi" << '\t' << "tau" << '\t' << "dq/dA" << '\t' << "dPhi/dA" << '\t' << "dtau/dA" << '\t';
    myfile << "d2q/dA2" << '\t' << "d2Phi/dA2" << '\t' << "d2tau/dA2" << '\t' << "d3tau/dA3" << std::endl;
    if(As.size() == 0) solve();
    for(int i = 0; i < As.size(); i++)
    {
        myfile << As[i] << '\t' << qs[i] << '\t' << Phis[i] << '\t' << taus[i] << '\t' << dqs[i] << '\t' << dPhis[i] << '\t' << dtaus[i] << '\t';
        myfile << d2qs[i] << '\t' << d2Phis[i] << '\t' << d2taus[i] << '\t' << d3taus[i] << std::endl; 
    }
    myfile.close();
}

void HVQCD::savePotentials(std::string path)
{
    /*
        Saves the values of the potentials Vg, Vf, k and w in as a function of A in the file
        specified by path.
        The data organized in columns in the following way: A  Vg  Vf  k  w
    */
    std::ofstream myfile;
    myfile.open(path);
    myfile << "A\tVg\tVf\tk\tw" << std::endl;
    // If the background has not been solved, solve it first
    if(As.size() == 0) solve();
    // Write the data into the file
    for(int i = 0; i < As.size(); i++)
    {
        myfile << As[i] << '\t' << Vg(Phis[i]) << '\t' <<  Vf(Phis[i], taus[i]) << '\t' << k(Phis[i]) << '\t' << w(Phis[i]) << std::endl;
    }
    myfile.close();
}

void HVQCD::saveZandUprofiles(std::string path)
{
    // Saves z and u profile in the file located at path
    std::ofstream myfile;
    myfile.open(path);
    myfile << "A" << '\t' << "z" << '\t' << "u" << std::endl;
    if(As.size() == 0) solve();
    for(int i = 0; i < As.size(); i++)
    {
        myfile << As[i] << '\t' << zs[i] << '\t' << us[i] << std::endl;
    }
    myfile.close();
}

void HVQCD::saveSchrodingerPotentials(std::string path)
{
    /* Computes and saves the Schrodinger potentials of Flavour Non-Singlet
       Vector, Axial vector, Pseudoscalar and Scalar mesons in a file.
       The text file is at path.
    */
    if(As.size() == 0) solve();
    // Compute the Schrodinger potentials
    std::vector<double> V2G = computeV2GPotential();
    std::vector<double> VVM = computeVectorMesonPotential();
    std::vector<double> VAVM = computeAxialVectorMesonNonSingletPotential(VVM);
    std::vector<double> VPSM = computePseudoScalarMesonPotential();
    std::vector<double> VSM = computeScalarMesonPotential();
    std::vector<double> VSingletAVM = computeAxialVectorMesonSingletPotential(VAVM);

    // Save the potentials in path
    std::ofstream myfile;
    myfile.open(path);
    myfile << "A" << '\t' << "2G" << '\t' << "VM" << '\t' << "AVM" << '\t' << "PSM" << '\t' << "SM" << '\t' << "SingAVM" << std::endl;
    for(int i = 0; i < As.size(); i++)
    {
        myfile << As[i] << '\t' << V2G[i] << '\t' << VVM[i] << '\t' << VAVM[i] << '\t' << VPSM[i] << '\t' << VSM[i] << '\t' << VSingletAVM[i] << std::endl;
    }
    myfile.close();

}

void HVQCD::showRatioValues(std::string method)
{
    /*
        This function starts by checking that the profile of the background fields has been computed.
        If not it will first solve the EOMs to determine the background fields.
        After the background fields are known it will compute and display the following ratios of the particle
        spectrum with respect to the mass of the lowest pseudoscalar (i.e. the pion)
    */
   // Check if the background fields have been computed
   if (As.size() == 0) this->solve();
   // Compute the Schrodinger potential of the fluctuations
    std::vector<double> V2G, VVM, VAVM, VPSM, VSM, VSingletAVM;
    if(add_tensor_glueball) V2G = computeV2GPotential();
    VVM = computeVectorMesonPotential();
    VAVM = computeAxialVectorMesonNonSingletPotential(VVM);
    if(add_scalars)
    {
        VPSM = computePseudoScalarMesonPotential();
        VSM = computeScalarMesonPotential();
    }
    if(add_singlet_axial) VSingletAVM = computeAxialVectorMesonSingletPotential(VAVM);
    // Compute the masses
    std::vector<double> TGMasses, VMMasses, AVMMasses, PSMMasses, SMMasses, SingletAVMMasses;
    if(add_tensor_glueball) TGMasses = computeMasses(zs, V2G, 1, method);
    VMMasses = computeMasses(us, VVM, 6, method);
    AVMMasses = computeMasses(us, VAVM, 5, method);
    if(add_scalars)
    {
        PSMMasses = computeMasses(us, VPSM, 5, method);
        SMMasses = computeMasses(us,VSM, 2, method);
    }
    if(add_singlet_axial) SingletAVMMasses = computeMasses(us, VSingletAVM, 2, method);
    // Compute the predicted ratios
    std::vector<double> RTGPred, RrhoPred, Ra1Pred, RpiPred, Ra0Pred, RomegaPred, Rf1Pred;
    if(add_tensor_glueball) RTGPred = {TGMasses[0]/ VMMasses[0]};
    RrhoPred = {VMMasses[1]/VMMasses[0], VMMasses[2]/VMMasses[0], VMMasses[3]/VMMasses[0], VMMasses[4]/VMMasses[0], VMMasses[5]/VMMasses[0]};
    Ra1Pred = { AVMMasses[0]/VMMasses[0], AVMMasses[1]/VMMasses[0], AVMMasses[2]/VMMasses[0], AVMMasses[3]/VMMasses[0], AVMMasses[4]/VMMasses[0]};
    if(add_scalars)
    {
        RpiPred = {PSMMasses[0]/VMMasses[0], PSMMasses[1]/VMMasses[0], PSMMasses[2]/VMMasses[0], PSMMasses[3]/VMMasses[0], PSMMasses[4]/VMMasses[0]};
        Ra0Pred = {SMMasses[0] / VMMasses[0], SMMasses[1] / VMMasses[0]};
    }
    if(add_singlet_vector) RomegaPred = {VMMasses[0]/VMMasses[0], VMMasses[1]/VMMasses[0], VMMasses[2]/VMMasses[0]};
    if(add_singlet_axial) Rf1Pred = {SingletAVMMasses[0]/VMMasses[0], SingletAVMMasses[1]/VMMasses[0]};
    // Compare the predicted ratios with the known ones
    std::cout << "Predicted Ratios" << '\t' << "Measured Ratios" << '\t' << "(Rpred - Robs) / Robs" << std::endl;
    // Tensor glueball ratio
     std::cout << "TENSOR GLUEBALL SECTOR" << std::endl;
    if(add_tensor_glueball) std::cout << RTGPred[0] << '\t' << RTG_rho[0] << '\t' << (RTGPred[0] - RTG_rho[0]) / RTG_rho[0] << std::endl;
    std::cout << std::endl;
    // Vector meson ratios
    std::cout << "VECTOR MESON NONSINGLET SECTOR" << std::endl;
    for(int i = 0; i < Rrho_rho.size(); i++) std::cout << RrhoPred[i] << '\t' << Rrho_rho[i] << '\t' << (RrhoPred[i] - Rrho_rho[i]) / Rrho_rho[i] << std::endl;
    std::cout << std::endl;
    // Axial vector meson ratios
    std::cout << "AXIAL VECTOR MESON NONSINGLET SECTOR" << std::endl;
    for(int i = 0; i < Ra1_rho.size(); i++) std::cout << Ra1Pred[i] << '\t' << Ra1_rho[i] << '\t' << (Ra1Pred[i] - Ra1_rho[i]) / Ra1_rho[i] << std::endl;
    std::cout << std::endl;
    if(add_scalars)
    {
        std::cout << "PSEUDOSCALAR MESON NONSINGLET SECTOR" << std::endl;
        for(int i = 0; i < Rpi_rho.size(); i++) std::cout << RpiPred[i] << '\t' << Rpi_rho[i] << '\t' << (RpiPred[i] - Rpi_rho[i]) / Rpi_rho[i] << std::endl;
        std::cout << std::endl;
        std::cout << "SCALAR MESON NONSINGLET SECTOR" << std::endl;
        for(int i = 0; i < Ra0_rho.size(); i++) std::cout << Ra0Pred[i] << '\t' << Ra0_rho[i] << '\t' << (Ra0Pred[i] - Ra0_rho[i]) / Ra0_rho[i] << std::endl;
        std::cout << std::endl;
    }
    // Singlet Vector meson ratios
    if(add_singlet_vector)
    {
        std::cout << "VECTOR MESON SINGLET SECTOR" << std::endl;
        for(int i = 0; i < Romega_rho.size(); i++) std::cout << RomegaPred[i] << '\t' << Romega_rho[i] << '\t' << (RomegaPred[i] - Romega_rho[i]) / Romega_rho[i] << std::endl;
        std::cout << std::endl;
    }
    // Singlet Axial vector meson ratios
    if(add_singlet_axial)
    {
        std::cout << "AXIAL VECTOR MESON SINGLET SECTOR" << std::endl;
        for(int i = 0; i < Rf1_rho.size(); i++) std::cout << Rf1Pred[i] << '\t' << Rf1_rho[i] << '\t' << (Rf1Pred[i] - Rf1_rho[i]) / Rf1_rho[i] << std::endl;
        std::cout << std::endl;
    }
    printQuarkMass();
}