#ifndef U1NNMODE_H
#define U1NNMODE_H

#include <vector>
#include <string>
#include "methods/interpolation/Poly_Interp.hpp"

void setupU1NNcomputation();

class U1NNMode
{
    private:
        double q2;                                           // Photon virtuality
        double * ISPACE, * FSPACE;                           // Pointers to arrays that contain all the info to compute fQ, dfQ and factor
        void copy(const U1NNMode &rhs);                      // copy method
        // Functions needed to compute the nonnormalizable modes
        static void f(double * X, double * Z, double *F, double * PARS);
        static void df(double * X, double * Z, double * DF, double * PARS);
        static void g(int * I, double * Z, double * G, double * PARS);
        static void dg(int * I, double * Z, double * DG, double * PARS);
        static void guess(double * X, double * Z, double * DMVAL);
        // Static auxiliary variables to compute the mode profile
        static std::vector<double> z;
        static Poly_Interp<double> t0;
        static Poly_Interp<double> t1;
        static const int IFDIM = 5000;
        static const  int IIDIM = 5000;
        friend void setupU1NNcomputation();
    public:
        U1NNMode(const double q2 = 1.0);                     // Class constructor
        U1NNMode(const U1NNMode &mode);                      // Copy constructor
        double Q2() const;
        double fQ(const double x) const;
        double dfQ(const double x) const ;
        double factor(const double x) const ;
        double Gsquared(const double x) const;
        void computeMode();                                 // Compute the non-normalizable modex
        void saveMode(std::string file_path = "") const;    // Crates a txt file with columns x, fQ, dfQ and fQ^2 + (dfQ/dx)^2
        U1NNMode& operator= (const U1NNMode &mode);         // Assignment operator
        bool operator< (const U1NNMode &mode) const;        // Operator < 
        bool operator> (const U1NNMode &mode) const;        // Operator >
        ~U1NNMode();                                        // Class destructor
};

#endif