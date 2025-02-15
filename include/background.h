#ifndef BACKGROUND_H
#define BACKGROUND_H

#include <vector>
#include <boost/numeric/ublas/vector.hpp>

 // state definition
typedef boost::numeric::ublas::vector< double > state;

class Background
{
    protected:
        double sc, VgIR;
        std::vector<double> qs, Phis, dqs, dPhis, d2Phis, As, zs;
        static const double lambda0;
        static const double V1;
        static const double V2;
        // Dilaton potential as a function of lambda
        double Vgl(const double l);
        // dVg/dPhi
        double dVgldlambda(const double l);
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
        // virtual copy function
        virtual void copy(const Background &bck);
        // operator() relevant for the EOMs
        virtual void eom(const state &X, state &dXdA, const double A);
        // operator() relevant for the observer
        virtual void observer(const state &X, const double A );
        virtual void solveRaw(const double AIR, const double AUV);
    protected:
        virtual void finalizeBackground(const double AIR, const double AUV);
    public:
        // Class constructor
        Background(const double ssc = 3.0, const double VVgIR = 2.05);
        // Class copy constructor
        Background(const Background & bck);
        // Class assignment opertor
        Background& operator=(const Background &rhs);
        virtual ~Background();
        void solve(const double AIR = -150, const double AUV = 50);
        // Getters of sc, VgIR, qs, Phis, dqs, dPhis, d2Phis, As and zs
        double get_sc() const;
        double get_VgIR() const;
        // Dilaton potential as a function of phi
        double Vg(const double phi);
        std::vector<double> q() const;
        std::vector<double> Phi() const;
        std::vector<double> dq() const;
        std::vector<double> dPhi() const;
        std::vector<double> d2Phi() const;
        std::vector<double> A() const;
        std::vector<double> z() const;
};

// Function that computes the potential of 2^{++} glueballs
std::vector<double> computeV2GPotential(const Background &bck);

std::vector<double> computeMasses(const std::vector<double> &z, const std::vector<double> &V, int nMasses, std::string method = "cheb");
#endif