#include <iostream>
#include "background.h"
#include "schrodinger/schrodinger.h"

// Definition of Background constructor
Background::Background(const double ssc, const double VVgIR): sc(ssc), VgIR(VVgIR){}

// Definition of Background copy constructor
Background::Background(const Background &bck):
    sc(bck.sc), VgIR(bck.VgIR), qs(bck.qs), Phis(bck.Phis),
    dqs(bck.dqs), dPhis(bck.dPhis), d2Phis(bck.d2Phis),
    As(bck.As), zs(bck.zs) {}


void Background::copy(const Background &rhs)
{
    // Copy all the data
    sc = rhs.sc;
    VgIR = rhs.VgIR;
    qs = rhs.qs;
    Phis = rhs.Phis;
    dqs = rhs.dqs;
    dPhis = rhs.dPhis;
    d2Phis = rhs.d2Phis;
    As = rhs.As;
    zs = rhs.zs;
}

// Definition of the assignment operator of Background class
Background& Background::operator=(const Background &rhs)
{   
    if (this == &rhs) return *this;
    // Copy all the data
    copy(rhs);
    return *this ;
}

// Definition of Background class destructor
Background::~Background() {}

const double Background::lambda0 = 8 * M_PI * M_PI;
const double Background::V1 = 44. / (9. * M_PI * M_PI);
const double Background::V2 = 4619. / (3888. * std::pow(M_PI,4.0));

double Background::Vgl(const double l)
{
    // Returns Vg as a function of lambda
    double ans = 12 + V1 * l + V2 * l * l / (1 + sc * l / lambda0 ) ;
    ans += 3 * std::exp(-lambda0/(sc * l)) * VgIR * std::pow(sc * l, 4.0/3) * std::sqrt(std::log(1+ sc*l/lambda0)) / (4 * std::pow(M_PI,8.0/3));
    return ans ;
}

double Background::dVgldlambda(const double l)
{
    // Returns dVg/dl as a function of lambda
    double ans = V1 + 2 * V2 * l / (1+sc*l/lambda0) - sc * V2 * l * l / (std::pow(1 + sc * l / lambda0,2.0) * lambda0);
    ans += 6.0 * std::exp(- lambda0 / (sc * l)) * sc * sc * VgIR * l * std::pow(sc * l / lambda0, 1.0/3) / ((1+sc*l/lambda0) * lambda0 * lambda0 * std::sqrt(std::log(1 + sc * l / lambda0)));
    ans += 12 * std::exp(- lambda0 / (sc * l)) * VgIR * std::pow(sc * l / lambda0, 1.0/3) * std::sqrt(std::log(1 + sc * l / lambda0)) / l ;
    ans += 12 * (4.0 / 3) * std::exp(- lambda0 / (sc * l)) * sc * VgIR * std::pow(sc * l / lambda0, 1.0/3) * std::sqrt(std::log(1 + sc * l / lambda0)) / lambda0;
    return ans;
}

double Background::Vg(const double phi)
{
    // Retuns the value of the dilaton potential given the value of the dilaton
    double l = std::exp(phi);
    return Vgl(l);
}

double Background::dVgdPhi(const double phi)
{
    // Give the value of the dilaton Phi returns dVg/dPhi
    // Compute dVg/dlambda, lambda = std::exp(phi)
    double l = std::exp(phi);
    double ans = dVgldlambda(l);
    // dVg/dPhi = lambda dVg/dlambda
    return l * ans;
}

double Background::d2VgdPhi2(double phi)
{
    // Returns d2Vg/d2Phi
    double l = std::exp(phi);
    double argl = sc*l/lambda0;
    double exparg = std::exp(-1/argl);
    double ans = V1*l + 4619*l*l/(972*std::pow(M_PI,4)*(1+argl)) + 4619*sc*sc*std::pow(l,4)/(1944*std::pow(M_PI,4)*std::pow(1+argl,3)*std::pow(lambda0,2)) - 23095*sc*std::pow(l,3)/(3888*std::pow(M_PI,4)*std::pow(1+argl,2)*lambda0)-3*exparg*sc*sc*VgIR*l*l*std::pow(argl,4./3)/(std::pow(1+argl,2)*std::pow(lambda0,2)*std::pow(std::log(1+argl),1.5));
    ans += 6*exparg*VgIR*std::pow(argl,4./3)/((1+argl)*std::sqrt(std::log(1+argl))) + 16*exparg*sc*sc*VgIR*l*l*std::pow(argl,1./3)/((1+argl)*std::pow(lambda0,2)*std::sqrt(std::log(1+argl))) - 6*exparg*sc*sc*VgIR*l*l*std::pow(argl,4./3)/(std::pow(1+argl,2)*std::pow(lambda0,2)*std::sqrt(std::log(1+argl)));
    ans += 6*exparg*sc*VgIR*l*std::pow(argl,4./3)*(1+1/argl)/((1+argl)*lambda0*std::sqrt(std::log(1+argl))) + 16*exparg*VgIR*std::pow(argl,1./3)*std::sqrt(std::log(1+argl)) + 16*exparg*sc*sc*VgIR*l*l*std::sqrt(std::log(1+argl))/(3*std::pow(argl,2./3)*std::pow(lambda0,2));
    ans += 12*exparg*VgIR*std::pow(argl,4./3)*lambda0*(-1+1/argl)*std::sqrt(std::log(1+argl))/(sc*l) + 16*exparg*sc*VgIR*l*std::pow(argl,1./3)*(1+1/argl)*std::sqrt(std::log(1+argl))/lambda0;
    return ans;
}

double Background::PhiIR(const double z)
{
    // Given z returns the dilaton value in the IR
    return -(23./16) - (151./(2304 * z * z)) + 1.5 * z * z - std::log(sc/lambda0) ;
}

double Background::AIR(const double z)
{
    // Given z returns the warp factor value A in the IR
    return (23.0/24) - (173 / (3456 * z * z)) - z * z + 0.25 * std::log(6 * z * z) - std::log(VgIR)/2.0;
}

double Background::dAIR(const double z)
{
    // Given z returns the value of dA/dz in the IR
    return 173 / (1728.0 * std::pow(z, 3.0)) + 1 / (2.0 * z) - 2 * z ;
}

double Background::dzdA(const double q, const double A)
{
    // Returns dz/dA = q(A) * exp(-A) given q and A
    return q * std::exp(-A);
}

double Background::dqYM(const double q, const double phi)
{
    // Returns dq/dA in the Yang-Mills theory
   return (12 * q - std::pow(q,3.0) * Vg(phi)) / 3.0 ;
}

double Background::d2qYM(const double q, const double phi)
{
    // Returns d2q/dA2 in the Yang-Mills theory 
   double vg = Vg(phi);
   double dvg = dVgdPhi(phi);
   double d2q = 16*q - (16./3)*std::pow(q,3)*vg+std::pow(q,5)*std::pow(vg,2)/3 +std::pow(q,3)*std::sqrt(1-std::pow(q,2)*vg/12)*dvg;
   return d2q;
}

double Background::dPhiYM(const double q, const double phi)
{
    // Returns dPhi/dA in the Yang-Mills theory
   return -0.5 * std::sqrt(3) * std::sqrt(12 - q * q * Vg(phi)) ;
}

double Background::d2PhiYM(const double q, const double phi)
{
    // Returns d2Phi/dA2 in the Yang-Mills theory
   double vg = Vg(phi);
   double dvg = dVgdPhi(phi);
   double ans = q * q * vg * std::sqrt(36.0 - 3.0 * q * q * vg)/6.0 - (3./8) * q * q * dvg;
   return ans;
}

void Background::eom(const state &X, state &dXdA, const double A)
{
    std::cout << "eom member function of Background class." << std::endl;
}

void Background::observer(const state &X, const double A)
{
    std::cout << "observer member function of Background class." << std::endl;
}

void Background::finalizeBackground(const double AIR, const double AUV)
{
    std::cout << "Background::finalizeBackground. Doing nothing" << std::endl;
}

void Background::solveRaw(const double AIR, const double AUV)
{
    std::cout << "Calling Background::solveRaw. Doing nothing." << std::endl; 
}

void Background::solve(const double AIR, const double AUV)
{
    solveRaw(AIR, AUV);
}

double Background::get_sc() const {return this->sc;}

double Background::get_VgIR() const {return this->VgIR;}

std::vector<double> Background::q() const {return this->qs;}

std::vector<double> Background::Phi() const {return this->Phis;}

std::vector<double> Background::dq() const {return this->dqs;}

std::vector<double> Background::dPhi() const {return this->dPhis;}

std::vector<double> Background::d2Phi() const {return this->d2Phis;}

std::vector<double> Background::A() const {return this->As;}

std::vector<double> Background::z() const {return this->zs;}

std::vector<double> computeV2GPotential(const Background &bck)
{
    // Compute the Schrodinger potential of the 2^{++} glueballs
    // Get A, q and dq values
    std::vector<double> As = bck.A(), qs = bck.q(), dqs = bck.dq();
    if (As.size() == 0) throw(std::runtime_error("Solve the holographic QCD vacuum first!"));
    // Compute V2G
    std::vector<double> V2G(As.size());
    for(int i = 0; i < As.size(); i++)
    {
        V2G[i] = (15.0/4) * std::exp(2*As[i]) / std::pow(qs[i],2.0) - 1.5 * std::exp(2 * As[i]) * dqs[i] / std::pow(qs[i], 3.0);
    }
    return V2G;
}

std::vector<double> computeMasses(const std::vector<double> &z, const std::vector<double> &V, int nMasses, std::string method)
{
    // Compute the spectrum
    List spectrum = computeSpectrum(z, V, 2 * nMasses, method);
    // The eigenvalues are mass squared so we need to take its square root
    std::vector<double> m2 = spectrum.Es;
    std::vector<double> m;
    int count = 0, i = 0;
    while ((count < nMasses) && (i < 2 * nMasses))
    {
        if (m2[i] >= 0)
        {
            m.push_back(std::sqrt(m2[i]));
            count++;
            i++;
        }
        else i++;
    }
    // Return masses
    return m;
}