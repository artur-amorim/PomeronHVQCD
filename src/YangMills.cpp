#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include "YangMills.h"
#include "schrodinger/common.h"
#include "methods/rootFind.hpp"
#include "methods/optimization/NelderMead.hpp"


// Definition YangMills constructor
YangMills::YangMills(const double ssc, const double VVgIR): Background(ssc, VVgIR) {}

// Definition YangMills copy constructor
YangMills::YangMills(const YangMills &ym): Background(ym) {}

// Definition YangMills assignment operator
YangMills& YangMills::operator=(const YangMills &rhs)
{   
    if (this == &rhs) return *this;
    // Copy all the data
    copy(rhs);
    d3Phis = rhs.d3Phis;
    return *this ;
}

// Definition YangMills destructor
YangMills::~YangMills(){}

double YangMills::d3PhiYM(const double q, const double phi)
{
    // Returns d3Phi/dA3 in the Yang-Mills theory
   double vg = Vg(phi);
   double dvg = dVgdPhi(phi);
   double d2vg = d2VgdPhi2(phi);
   double sqrtarg = std::sqrt(1-q*q*vg/12.) ;
   double ans = 8*q*q*vg/sqrtarg - (5./3)*std::pow(q*q*vg,2)/sqrtarg + std::pow(q*q*vg,3)/(12*sqrtarg);
   ans += -6*q*q*dvg + (5./8)*std::pow(q,4)*vg*dvg + (9./8)*q*q*d2vg/sqrtarg - (3./32)*std::pow(q,4)*vg*d2vg/sqrtarg;
   return ans;
}

// Definition of Yang-Mills equations of motion
void YangMills::eom(const state &X, state &dXdA, const double A)
{
    // X = q, Phi, z
    dXdA[0] = dqYM(X[0], X[1]);
    dXdA[1] = dPhiYM(X[0], X[1]);
    dXdA[2] = dzdA(X[0], A);
}

// Definition of observer function
void YangMills::observer(const state &X, const double A)
{
    // X = (q, Phi, z)
    double dq = dqYM(X[0], X[1]);
    double dphi = dPhiYM(X[0], X[1]);
    double d2q = d2qYM(X[0], X[1]);
    double d2phi = d2PhiYM(X[0], X[1]);
    qs.push_back(X[0]);
    Phis.push_back(X[1]);
    dqs.push_back(dq);
    dPhis.push_back(dphi);
    d2Phis.push_back(d2phi);
    d3Phis.push_back(d3PhiYM(X[0], X[1]));
    As.push_back(A);
    zs.push_back(X[2]);
}

void YangMills::finalizeBackground()
{
    // Selects the relevant values to compute the potentials later
    std::vector<double> A, z, q, Phi, dq, dPhi, d2Phi, d3Phi;
    for(int i = 0; i < As.size(); i++)
    {
        // For the spectrum we are only interested in -50 < A < 20
        if ((zs[i]-zs.back()> 1e-6) && (Phis[i] < 120))
        {
            A.push_back(As[i]); z.push_back(zs[i]-zs.back());
            q.push_back(qs[i]); Phi.push_back(Phis[i]);
            dq.push_back(dqs[i]); dPhi.push_back(dPhis[i]);
            d2Phi.push_back(d2Phis[i]); d3Phi.push_back(d3Phis[i]);
        }
    }
    As = A; zs = z; qs = q; Phis = Phi, dqs = dq; dPhis = dPhi;
    d2Phis = d2Phi; d3Phis = d3Phi;
}


// Solve Yang-Mills equations of motion
void YangMills::solve()
{
    // Setting up IR boundary conditions
   double Air = -150, AUVYM = 50, h = 0.01;
   std::function<double(double)> func = [this, &Air] (double z) { return this->AIR(z) - Air ;} ;
   double zIRYM = zbrent(func, 0.1, 100.0, 1e-9, true);
   double qir = exp(Air) / dAIR(zIRYM), Phiir = PhiIR(zIRYM);
   state X(3);
   X <<= qir, Phiir, zIRYM;
   // Solve the EOMs
   auto eomfun = [this] (const state &Y, state &dYdA, double A) {this->eom(Y, dYdA, A);};
   auto obsfun = [this] (const state &Y, double A) {this->observer(Y, A);};
   typedef boost::numeric::odeint::result_of::make_dense_output<boost::numeric::odeint::runge_kutta_dopri5<state> >::type dense_stepper;
   dense_stepper stepper = make_dense_output(1.0e-12, 1.0e-12, boost::numeric::odeint::runge_kutta_dopri5< state >());
   integrate_const(stepper, eomfun, X, Air, AUVYM, h, obsfun);
   finalizeBackground();
}

std::vector<double> YangMills::d3Phi() const {return this->d3Phis ;}

std::vector<double> computeV0GPotential(const YangMills &ym)
{
    std::vector<double> As = ym.A(), qs = ym.q(), dqs = ym.dq(), dPhis = ym.dPhi(), d2Phis = ym.d2Phi(), d3Phis = ym.d3Phi();
    if (As.size() == 0) throw(std::runtime_error("Solve the holographic QCD vacuum first!"));
    std::vector<double> V0G(As.size());
    for(int i = 0; i < As.size(); i++)
    {
        V0G[i] = 3.75 - 1.5 * dqs[i] / qs[i] + 4 * d2Phis[i]/dPhis[i] - dqs[i] * d2Phis[i] / (qs[i] * dPhis[i]) + d3Phis[i] / dPhis[i] ;
        V0G[i] = std::exp(2*As[i]) * V0G[i] / std::pow(qs[i],2) ;
    }
    return V0G;
}

// Computes the Yang-Mills theory spectrum
void computeYangMillsSpectrum(const YangMills &ym, const int n_tensor, const int n_scalar)
{
    std::vector<double> zs = ym.z();
    // Compute tensor glueball spectrum
    std::vector<double> VG2 = computeV2GPotential(ym);
    std::vector<double> tensor_masses = computeMasses(zs, VG2, n_tensor, "cheb");
    // Compute scalar glueball spectrum
    std::vector<double> VG0 = computeV0GPotential(ym);
    std::vector<double> scalar_masses = computeMasses(zs, VG0, n_scalar, "cheb");

    // Display the tensor glueball masses
    std::cout << "TENSOR GLUEBALL MASSES - YANG-MILLS" << std::endl;
    for(int i = 0; i < n_tensor; i++) std::cout << tensor_masses[i] << '\t';
    std::cout << std::endl;
    std::cout << std::endl;

    // Display the scalar glueball masses
    std::cout << "SCALAR GLUEBALL MASSES - YANG-MILLS" << std::endl;
    for(int i = 0; i < n_scalar; i++) std::cout << scalar_masses[i] << '\t';
    std::cout << std::endl;
    std::cout << std::endl;
}

// Computes the Yang-Mills theory glueball ratios
void computeYangMillsRatios(const YangMills &ym)
{
    // Displays the ratios R20 = m_{2++} / m_{0++} and R00 = m_{0*++} / m_{0++}

    std::vector<double> zs = ym.z();
    // Compute tensor glueball spectrum
    std::vector<double> VG2 = computeV2GPotential(ym);
    std::vector<double> tensor_masses = computeMasses(zs, VG2, 1, "cheb");
    // Compute scalar glueball spectrum
    std::vector<double> VG0 = computeV0GPotential(ym);
    std::vector<double> scalar_masses = computeMasses(zs, VG0, 2, "cheb");

    std::cout << "YANG MILLS GLUEBALL RATIOS" << std::endl;
    // Display the R20 ratio
    std::cout << "R20: " << tensor_masses[0] / scalar_masses[0] << std::endl;
    // Display the R00 ratio
    std::cout << "R00: " << scalar_masses[1] / scalar_masses[0] << std::endl;
    std::cout << std::endl;
    
}

// Fits the parameters to the Yang-Mills ratios
void fitYangMills(const double sc_guess, const double VgIR_guess)
{
    // Define function to be fitted
    auto J = [] (const std::vector<double> X)
    {
        double sc = X[0], VgIR = X[1];
        // Create Yang-Mills object
        YangMills ym(sc, VgIR);
        // Solve Yang-Mills theory
        ym.solve();

        // Find the glueball spectrum
        std::vector<double> zs = ym.z();
        // Compute tensor glueball spectrum
        std::vector<double> VG2 = computeV2GPotential(ym);
        std::vector<double> tensor_masses = computeMasses(zs, VG2, 1, "cheb");
        // Compute scalar glueball spectrum
        std::vector<double> VG0 = computeV0GPotential(ym);
        std::vector<double> scalar_masses = computeMasses(zs, VG0, 2, "cheb");

        double ans = 0;

        // Add relative difference squared of R00
        ans += std::pow((scalar_masses[1]/scalar_masses[0] - R00)/R00,2);
        // Add relative difference squared of R20
        ans += std::pow((tensor_masses[0]/scalar_masses[0] - R20)/R20,2);
        std::cout << "sc: " << sc << '\t' << "VgIR: " << VgIR << '\t' << "J: " << ans << std::endl;
        return ans;

    };

    // Define guess point
    std::vector<double> X_guess = {sc_guess, VgIR_guess};
    // Fit the model
    std::vector<double> X_optimum = optimFunction(X_guess, J, 0.5);

    // Display the results
    std::cout << "Values of sc and VgIR for the best fit" << std::endl;
    std::cout << "sc: " << X_optimum[0] << '\t' << "VgIR: " << X_optimum[1] << std::endl;

    // Create Yang-Mills object and display the ratios
    YangMills ym(X_optimum[0], X_optimum[1]);
    ym.solve();
    computeYangMillsRatios(ym);
}
