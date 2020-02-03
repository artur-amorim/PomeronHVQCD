#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include "YangMills.h"
#include "schrodinger/common.h"


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
    return *this ;
}

// Definition YangMills destructor
YangMills::~YangMills(){}

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
    As.push_back(A);
    zs.push_back(X[2]);
}

void YangMills::finalizeBackground()
{
    // Make zUV = 0
    double zUV = zs.back();
    for(int i = 0; i < zs.size(); i++) zs[i] = zs[i] - zUV;
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