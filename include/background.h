#ifndef BACKGROUND_H
#define BACKGROUND_H

#include <vector>
#include <boost/numeric/ublas/vector.hpp>

class Background
{
    protected:
        // state definition
        typedef boost::numeric::ublas::vector< double > state;
        double sc, VgIR;
        std::vector<double> qs, Phis, dqs, dPhis, d2Phis, As, zs;
        static const double lambda0;
        static const double V1;
        static const double V2;
        // Dilaton potential as a function of lambda
        double Vgl(const double l);
        // dVg/dPhi
        double dVgldlambda(const double l);
        // Dilaton potential as a function of phi
        double Vg(const double phi);
        // dVg/dPhi
        double dVgdPhi(const double phi);
        // d2Vg/dPhi2
        double d2VgdPhi2(const double phi);
        // Dilaton profile in the IR
        double PhiIR(const double z);
        // Warp factor profile in the IR
        double AIR(const double z);
        double dAIR(const double z);
        double dzdA(const double q, const double A);
        // Relevant derivatives of the background fields
        double dqYM(const double q, const double phi);
        double d2qYM(const double q, const double phi);
        double dPhiYM(const double q, const double phi);
        double d2PhiYM(const double q, const double phi);
        // operator() relevant for the EOMs
        virtual void eom(const state &X, state &dXdA, const double A);
        // operator() relevant for the observer
        virtual void observer(const state &X, const double A );
        virtual void finalizeBackground();
    public:
        Background(const double ssc = 3.0, const double VVgIR = 2.05);
        virtual void solve();
        virtual std::vector<std::vector<double> > getBackgroundFields();
        virtual std::vector< double > getBackgroundPars();
        std::vector<double> computeV2GPotential();
        // Declaration of computeSpectrum function of Background class
        virtual void computeSpectrum();
        virtual double J(const state &x);
};

std::vector<double> computeMasses(const std::vector<double> &z, const std::vector<double> &V, int nMasses, std::string method = "cheb");
#endif