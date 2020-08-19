#include <iostream>
#include <fstream>
#include <stdexcept>
#include <boost/algorithm/string.hpp>
#include "HolographicVQCD.h"
#include "Process.h"
#include "Sigma.h"
#include "Reggeon.h"
#include "Spectra.h"
#include "methods/interpolation/Poly_Interp.hpp"
#include "methods/vectorOperators.hpp"
#include "methods/search.hpp"

extern"C"
{
    void dqags_(double (*f)(double *, void *), void * params,
                double * a, double * b, double * epsabs, double * epsrel, 
                double * result, double * abserr, int * neval,int * ier,
                int * limit, int * lenw, int * last, int * iwork, double * work);
}

class SigmaIzNIntegrand
{
    private:
        Poly_Interp<double> f1, f2, f3;
    public:
        SigmaIzNIntegrand(const Poly_Interp<double> &func1, const Poly_Interp<double> &func2, const Poly_Interp<double> &func3);
        double operator()(const double x);
};

SigmaIzNIntegrand::SigmaIzNIntegrand(const Poly_Interp<double> &func1, const Poly_Interp<double> &func2, const Poly_Interp<double> &func3): 
                                    f1(func1), f2(func2), f3(func3) {}

double SigmaIzNIntegrand::operator()(const double x) {return f1.interp(x) * f2.interp(x) * f3.interp(x) ;}

double fSigma(double * x, void * params)
{
    return ((SigmaIzNIntegrand *) params)->operator()(*x);
}

void Sigma::loadData(const std::string & file_path, const double conv_factor)
{
    /*
        This function loads the data of sigma contained in the path file_path
        The data is assumed to be in 6 columns:
        W(GeV), WErrorp, WErrorm, sigma(nb, mub, mb), Errorp and Errorm
        WErrorp is the positive error uncertainty in W
        WErrorm is the negative error uncertainty in W
        sigma(nb, mub, mb) is the cross section in nb, mub or mb
        Errorp is the plus error in sigma
        Errorm is the minus error in sigma 
    */
   // Check that the file path is not empty
   if (file_path.size() == 0)
   {
       std::cout << "No file path given. Exiting SigmaGammaGamma::loadData" << std::endl;
       return;
   }
   std::ifstream file;
   file.open(file_path);
   if(file.fail()) std::runtime_error("Error: File with sigma(gamma p -> hadrons) data not opened."); 
   std::string line ;
   getline(file, line) ;
   std::cout << "Loading sigma data" << std::endl;
   std::vector<double> Ws, WsPlus, WsMinus, sigmas, sigmaErrs;
   std::vector<std::string> result;
   // Convert GeV to U units
   const double mrho_U_mrho_GeV = 4.3669;
   while(getline(file, line))
    {
        boost::split(result, line, boost::is_any_of("\t") ) ;
        Ws.push_back(stod(result[0]) * mrho_U_mrho_GeV) ;
        WsPlus.push_back((stod(result[0]) + stod(result[1])) * mrho_U_mrho_GeV);
        WsMinus.push_back((stod(result[0]) - stod(result[2])) * mrho_U_mrho_GeV);
        sigmas.push_back(stod(result[3]) * conv_factor);
        sigmaErrs.push_back(std::max(stod(result[4]), stod(result[5])) * conv_factor);
    }
    // Check that all std::vector containers have the same side
    if( Ws.size() != sigmas.size() || Ws.size() != sigmaErrs.size() || sigmas.size() != sigmaErrs.size())
    {
        throw std::runtime_error("Ws, sigma and sigmaErrs vector containers in SigmaGammaProton don't have the same size"); 
    }
    this->setDataPts({Ws, WsPlus, WsMinus, sigmas, sigmaErrs});
}

void Sigma::copy(const Sigma &sigma)
{
    // Auxiliary copy function of the Sigma class
    // Copy data points using the Process copy
    Process::copy(sigma);
    // Finally copy barn_to_UMINUS2
    barn_to_UMINUS2 = sigma.barn_to_UMINUS2;
}

// Define the class constructor
Sigma::Sigma(const std::string &file_path, const double conv_factor) : Process(), barn_to_UMINUS2(conv_factor)
{
    loadData(file_path, barn_to_UMINUS2);
}

// Define the class copy constructor
Sigma::Sigma(const Sigma &sigma): Process(sigma), barn_to_UMINUS2(sigma.barn_to_UMINUS2) {}

// Define expVal
std::vector<double> Sigma::expVal()
{
    // Returns a std::vector container with the values of sigma
    return this->getDataPts()[3];
}

// Define expErr
std::vector<double>  Sigma::expErr()
{
    // Returns a std::vector container with the values of sigmaErrs
    return this->getDataPts()[4];
}

std::vector<std::vector<double> >  Sigma::expKinematics()
{
    /*
        Returns a std::vector<std::vector<double> > with the values of Ws.
        This format is just to make it uniform with processed that depend
        in more than one kinematical variable
    */
   std::vector<double> Ws = this->getDataPts()[0], WsPlus = this->getDataPts()[1], WsMinus = this->getDataPts()[2];
   return {Ws, WsPlus, WsMinus};
}

double Sigma::IzN(const std::vector<double> &kin, const Reggeon &reg)
{
    /*
        Computes the IzN integral that appears in the 
        holographic computation of sigma(gamma proton -> hadrons)
        kin - std::vector<double> with just one element: W
        reg - Reggeon object from wich j_n and \psi_n can be accessed.
        returns \int du e^(-(j_n - 3/2)A) exp^(-7/3 Phi) * Vf * ws^2 \psi_n if the reggeon comes from the pomeron trajectory
        or returns \int du e^{-(j_n - 2) A_s} sqrt(exp^(-10\Phi/3) * Vf * ws^2) \psi_n if it comes from the meson trajectory

    */
    // Get J and the wavefunction wf from Reggeon object reg
    const double J = reg.getJ();
    const std::string reg_name = reg.getName();
    std::vector<std::vector<double> > wf = reg.getWf();

    // Define e^(As(1.5 - jn)) or e^(As(2-j_n))
    std::vector<double> fact1;
    if(reg_name == "gluon") fact1 = exp((1.5-J) * Astring);
    else fact1 = exp((2.0-J) * Astring);
    Poly_Interp<double> f1(u, fact1, 4);

    // Choose the background potentials factor accordingly
    Poly_Interp<double> bckPotFac;
    if(reg_name == "gluon") bckPotFac = potFactor;
    else bckPotFac = MesonPotFactor;

    // Interpolate the wavefunction after computing the respective u value.
    if(reg_name == "gluon")  for(int i = 0; i < wf[0].size(); i++) wf[0][i] = ufunc.interp(wf[0][i]);
    Poly_Interp<double> f3(wf[0], wf[1], 4);


    // Define the integrand object
    SigmaIzNIntegrand integrand(f1, bckPotFac, f3);
    void * params = &integrand;
    // Compute the integral
    double izn = 0.0, abserr = 0.0;
    // If we look at the potential profiles vs u we see that this is good enough
    double a = u[0], b = 4;
    // Absolute and relative tolerances desidered
    double epsabs = 1e-9, epsrel = 1e-9;
    // Output variables 
    int neval = 0, ier;
    // Setup the workspace for the method
    int limit = 10000;
    int lenw = 100000;
    int last = 0;
    int * iwork = new int[limit];
    double * work = new double[lenw];
    // Evaluate the integral
    dqags_(fSigma, params, &a, &b, &epsabs, &epsrel, &izn, &abserr, &neval, &ier, &limit, &lenw, &last, iwork, work);
    // Free workspace from memory
    delete[] iwork;
    delete[] work;
    return izn;
}

double Sigma::IzNBar(const std::vector<double> &kin, const Reggeon &reg, const std::vector<double> &gs)
{
    /*
        Computes the IzNBar integral that appears in the 
        holographic computation of sigma(gamma gamma -> hadrons)
        kin - std::vector<double> with just one element: W
        reg - Reggeon object from wich the reggeon index can be accessed.
        gs - std::vector<double> that contains the constant quantities in our problem.
        returns s^(j_n -1 ) g_n associated with reggeon n
    */
   // Compute s
   double s = std::pow(kin[0],2), J = reg.getJ();
   const int reg_index = reg.getIndex();
   double iznbar = std::pow(s, J - 1) * gs[reg_index-1];
   return iznbar;
}

std::vector<kinStruct>  Sigma::getIzs(const std::vector< std::vector<double> > &points, const std::vector<Spectra> &spec)
{
    /*
        Computes all the Izs relevant to sigma(gamma p -> hadrons)
    */
    // Get the Reggeons in spec
    const std::vector<Reggeon> reggeons = spec[0].getReggeons();
    const int n_reggeons = reggeons.size();
    // Create vector of kinStructs
    std::vector< kinStruct > ans(1);
    std::vector<double> izns(n_reggeons,0);
    for(int reg_idx = 0; reg_idx < n_reggeons; reg_idx++) izns[reg_idx] = IzN({}, reggeons[reg_idx]);
    kinStruct izs(0.0, izns);
    ans[0] = izs;
    return ans ;
}

std::vector<kinStruct> Sigma::getIzsBar(const std::vector< std::vector<double> > &points, const std::vector<Spectra> &spec,
                                          const std::vector<double>  &gs)
{
    /*
        Computes all the IzNbars relevant to sigma(gamma proton -> hadrons)
    */
    // Create list of unique Ws
    std::vector<double> Ws = points[0], WsPlus = points[1], WsMinus = points[2];
    const int n_Ws = Ws.size();
    std::vector<kinStruct> ans(n_Ws) ;
    // Get kernels in spec. For sigma observables there is only one value of t, i.e. 0
    std::vector<Reggeon> reggeons = spec[0].getReggeons();
    const int n_reggeons = reggeons.size();
    std::vector<double> iznbars(3 * n_reggeons,0);
    for(int i = 0; i < n_Ws; i++)
    {
        for(int reg_index = 0; reg_index < n_reggeons; reg_index++)
        {
            iznbars[reg_index] = IzNBar({Ws[i]}, reggeons[reg_index], gs);
            iznbars[n_reggeons + reg_index] = IzNBar({WsPlus[i]}, reggeons[reg_index], gs);
            iznbars[2*n_reggeons + reg_index] = IzNBar({WsMinus[i]}, reggeons[reg_index], gs);
        }
        kinStruct izbars(Ws[i], iznbars);
        ans[i] = izbars;
    }
    return ans ;
}

std::vector<double>  Sigma::predict(const std::vector<kinStruct>  &Izs, const std::vector<kinStruct>  &IzsBar,
                                     const std::vector< std::vector<double> >  &points, const bool savePredictions)
{
    // Check that Izs, IzsBar and points have the same size
    if (points.size() == 0) throw std::runtime_error("points has length 0. Aborting Sigma::predict.");
    std::vector<double> ans(points[0].size()) ;
    kinStruct iznStruct, iznbarStruct;
    std::vector<double> izn, iznbar, central_value_iznbar;
    for(int i = 0; i < points[0].size(); i++)
    {
        iznStruct = binary_search<kinStruct>(Izs, kinStruct(0.0,{}));
        iznbarStruct = binary_search<kinStruct>(IzsBar, kinStruct(points[0][i],{}));
        izn = iznStruct.izns; iznbar = iznbarStruct.izns;
        // IzNBar has three times the number of elements of izn.
        // To predict the cross section we only need the IzBars computed at W
        // To do that we select the first izn.size() elements
        central_value_iznbar = {};
        for(int j = 0; j < izn.size(); j++) central_value_iznbar.push_back(iznbar[j]);
        ans[i] = sum(izn * central_value_iznbar);
    }
    if(savePredictions)
    {
        const double mrho_U_mrho_GeV = 4.3669;
        std::ofstream myfile;
        std::string file_path;
        std::cout << "Please introduce the path to save the predictions of sigma" << std::endl;
        std::cin >> file_path;
        myfile.open(file_path);
        myfile << "W\tPred" << std::endl;
        for(int i = 0; i < points[0].size(); i++) myfile << points[0][i] / mrho_U_mrho_GeV << '\t' << ans[i] / barn_to_UMINUS2 << std::endl;
    }
    return ans ;
}

std::vector<double> Sigma::diffObsWeighted(const std::vector<kinStruct> &Izs, const std::vector<kinStruct> &IzsBar, const std::vector< std::vector < double > > &points)
{
    if( points.size() == 0) std::vector< std::vector< double > > points = this->expKinematics() ;   // If points is NULL provide the experimental ones
    std::vector<double> Opred(points[0].size()), OWPlus(points[0].size()), OWMinus(points[0].size()) ;
    kinStruct iznStruct, iznbarStruct;
    std::vector<double> izn, iznbar, W_iznbar, Wplus_iznbar, Wminus_iznbar;
    for(int i = 0; i < points[0].size(); i++)
    {
        iznStruct = binary_search<kinStruct>(Izs, kinStruct(0.0,{}));
        iznbarStruct = binary_search<kinStruct>(IzsBar, kinStruct(points[0][i],{}));
        izn = iznStruct.izns; iznbar = iznbarStruct.izns;
        // IzNBar has three times the number of elements of izn.
        // To predict the cross section we only need the IzBars computed at W
        // To do that we select the first izn.size() elements
        const int izn_size = izn.size();
        W_iznbar = {}; Wplus_iznbar = {}; Wminus_iznbar = {};
        for(int j = 0; j < izn_size; j++)
        {
            W_iznbar.push_back(iznbar[j]);
            Wplus_iznbar.push_back(iznbar[izn_size + j]);
            Wminus_iznbar.push_back(iznbar[2 * izn_size + j]);
        }
        Opred[i] = sum(izn * W_iznbar);
        OWPlus[i] = sum(izn * Wplus_iznbar);
        OWMinus[i] = sum(izn * Wminus_iznbar);
    }
    const std::vector<double> Oexp  = this->expVal() ;                                              // Experimental values of the process
    std::vector<double> Oerr  = this->expErr() ;                                                    // Experimental errors of the process
    // Because we have uncertainty in W we need to add the effective uncertainty
    const std::vector<double> Oeff_uncert = maximum(abs(Opred-OWPlus), abs(Opred - OWMinus) );
    Oerr = sqrt(Oerr * Oerr + Oeff_uncert * Oeff_uncert) ;
    return (Opred - Oexp) / Oerr ;

}


Sigma::~Sigma()
{
    // Class destructor
}