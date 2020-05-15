#include <fstream>
#include <stdexcept>
#include <boost/algorithm/string.hpp>
#include "HolographicVQCD.h"
#include "DeepInelasticScattering.h"
#include "methods/search.hpp"
#include "methods/vectorOperators.hpp"
#include "methods/interpolation/Spline_Interp.hpp"

std::vector<double> DeepInelasticScattering::z = {};
std::vector<double> DeepInelasticScattering::Astring = {};
std::vector<double> DeepInelasticScattering::Vfw2fac = {};
Poly_Interp<double> DeepInelasticScattering::potFactor({},{},4);


DeepInelasticScattering::DeepInelasticScattering(const bool rrsslog, std::string file_path): ProcessObservable(rrsslog), modes({})
{
    // Load data
    loadData(file_path);
    // Compute all the necessary U(1) NN modes
    computeU1NNModes();
    // We now compute potFactor
    if (z.size() == 0)
    {
        z = hvqcd().z();
        Astring = hvqcd().Astring(); 
        Vfw2fac.resize(z.size());
        // e^{-7/3 \Phi} V_f w_s^2 = e^(\Phi / 3) V_f w^2 in terms of the Einstein frame background potentials
        std::vector<double> Phis = hvqcd().Phi(), taus = hvqcd().tau();
        for(int i = 0; i < z.size(); i++) Vfw2fac[i] = std::exp(Phis[i] / 3) * hvqcd().Vf(Phis[i], taus[i]) * std::pow(hvqcd().w(Phis[i]),2) ;
        // Now we reverse z, Astring and Vfw2fac because we want them from the UV to the IR
        std::reverse(std::begin(z), std::end(z));
        std::reverse(std::begin(Astring), std::end(Astring));
        std::reverse(std::begin(Vfw2fac), std::end(Vfw2fac));
        potFactor = Poly_Interp<double>(z, Vfw2fac, 4);
    }
}

DeepInelasticScattering::DeepInelasticScattering(const DeepInelasticScattering &dis): ProcessObservable(dis), modes(dis.modes){}

void DeepInelasticScattering::loadData(std::string file_path)
{
     // Vector container with Q2, x, Fi and Fierr values, with i = 2, L
    if (file_path.size() == 0) return;
    std::ifstream file;
    file.open(file_path);
    if(file.fail()) std::runtime_error("Error: File with DIS data not opened."); 
    std::string line ;
    getline(file, line) ;
    std::cout << "Loading DIS data" << std::endl;
    std::vector<double> Q2, x, Fi, Fierr;
    std::vector<std::string> result;
    // Ratio of the rho mass in the model to the rho mass in GeV. Assuming we are using global fit results to HVQCD
    const double mrho_U_mrho_GeV = 4.06789;
    while(getline(file, line))
    {
        boost::split(result, line, boost::is_any_of("\t") ) ;
        Fi.push_back(std::stod(result[0])) ;
        Q2.push_back(std::stod(result[1]) * std::pow(mrho_U_mrho_GeV, 2) ) ;
        x.push_back(std::stod(result[2])) ;
        Fierr.push_back(stod(result[3]));
    }
    // Check that all std::vector containers have the same side
    if( Fi.size() != Q2.size() || Q2.size() != x.size() || x.size() != Fierr.size())
    {
        std::cout << "Vector containers don't have the same size" << std::endl;
        throw "error";
    }
    this->setDataPts({Q2, x, Fi, Fierr}) ;
}

void DeepInelasticScattering::computeU1NNModes()
{
    /*
        Computes the relevant U(1) Nonnormalizable modes from the know values of Q2
        present in the data points.
        The modes are stored in the std::vector container modes and are ordered by Q2.
    */
    // First clear the modes container
    modes.clear();
    // Get a std::vector with a unique list of Q2
    std::vector<std::vector<double> > data = this->getDataPts();
    std::vector<double> Q2s(0);
    if(data.size()!= 0) Q2s = data[0];
    else return ;
    std::sort(Q2s.begin(), Q2s.end());
    Q2s.erase(unique(Q2s.begin(), Q2s.end()), Q2s.end());
    // For each Q2 compute a mode
    std::cout << "Computing the U(1) Nonnormalizable modes of DIS object:" << std::endl;
    double Q2;
    for(int i = 0; i < Q2s.size(); i++)
    {
        Q2 = Q2s[i];
        std::cout << "Computing mode for Q2 = " << Q2 << std::endl;
        // Compute the mode
        U1NNMode mode(Q2);
        mode.computeMode();
        modes.push_back(mode);
    }
}

void DeepInelasticScattering::copy(const DeepInelasticScattering &rhs)
{
    ProcessObservable::copy(rhs);
    modes = rhs.modes;
}

U1NNMode DeepInelasticScattering::searchMode(const double Q2)
{
    /*
        Function that searches for the mode with virtuality Q2.
        Uses the binary search algorithm to find the mode.
        If the mode is not found it is computed.
    */
   U1NNMode mode(Q2);
   try
   {
       return binary_search<U1NNMode>(modes, mode);
   }
   catch(const char *)
   {
       std::cout << Q2 << std::endl;
       mode.computeMode();
       return mode;
   }
}

std::vector<U1NNMode> DeepInelasticScattering::getModes()
{
    /*
        Returns all the relevant U(1) nonnormalizable modes.
        If the modes have not been computed it computes the modes first.
    */

    if(modes.size()==0) this->computeU1NNModes();
    return this->modes;
}

std::vector<double> DeepInelasticScattering::expVal()
{
    // F2 data is present in index 2 of datapts
    return this->getDataPts()[2] ;
}


std::vector<double> DeepInelasticScattering::expErr()
{
    // F2 err is present in index 3 of datapts
    return this->getDataPts()[3] ;
}


std::vector < std::vector<double> > DeepInelasticScattering::expKinematics()
{
    std::vector < std::vector < double > > ans ;
    // Get the Q2s std::vector
    ans.push_back(this->getDataPts()[0]) ;
    // Get the xs std::vector
    ans.push_back(this->getDataPts()[1]) ;
    return ans ;
}

std::vector<double> DeepInelasticScattering::getNeededTVals()
{
    /*
        Returns a vector with a single zero entry because in the
        holographic computation of F2 and FL we only need t = 0
    */
    std::vector<double> ans(1, 0.0) ;
    return ans ;
}

double DeepInelasticScattering::IzNBar(std::vector<double> kin, const Reggeon &reg, const std::vector<double> &gs)
{
    /*
        Given a reggeon and it the kinematical point (Q2, x),
        it returns the right gn from the list gs and multiplies
        it by x^(1-J)
    */
   const double alpha0 = 0.0072973525693;
   const double x = kin[0];
   const double J = reg.getJ();
   const int reg_index = reg.getIndex();
   double iznbar = std::pow(x, 1-J) * gs[reg_index-1] / (4*M_PI*M_PI*alpha0);
   return iznbar;
}

std::vector<kinStruct> DeepInelasticScattering::getIzs(const std::vector< std::vector<double> > &points,
                                                                   const std::vector<Spectra> &spec)
{
    // Create unique list of Q2s
    std::vector<double> Q2s = points[0];
    sort(Q2s.begin(), Q2s.end());
    Q2s.erase(unique(Q2s.begin(), Q2s.end()), Q2s.end());
    const int n_Q2s = Q2s.size();
    // Get the Reggeons in spec
    const std::vector<Reggeon> reggeons = spec[0].getReggeons();
    const int n_reggeons = reggeons.size();
    // Compute the modes first if not computed previously
    if(modes.size() == 0) this->computeU1NNModes();
    // Create vector of kinStructs
    std::vector< kinStruct > ans(n_Q2s);
    std::vector<double> kinematics, izns(n_reggeons,0);
    for(int i = 0; i < n_Q2s; i++)
    {
        kinematics = {Q2s[i]};
        for(int reg_idx = 0; reg_idx < n_reggeons; reg_idx++) izns[reg_idx] = IzN(kinematics, reggeons[reg_idx]);
        kinStruct izs(Q2s[i], izns);
        ans[i] = izs;
    }
    return ans ;
}


std::vector<kinStruct> DeepInelasticScattering::getIzsBar(const std::vector< std::vector<double> > &points,
                                                                      const std::vector<Spectra>& spec,
                                                                      const std::vector<double> &gs)
{
    /*
        Computes all the necessary IzNBar.
        points: std::vector of vectors. Each std::vector contains a pair (Q2, x).
        spec: object which contains a std::vector with different kernels.
        gs: std::vector with the gn corresponding to each reggeon present in spec.
        To get the correct gs for each point we initiate the variable gs_index.
        Then it iterates over the kernels and then over the reggeons of each kernel.
        In each iteration gs_index is increased by one.
        This guarantees that the correct gn is obtained.        
    */
    // Create list of unique xs
    std::vector<double> xs = points[1];
    std::sort(xs.begin(), xs.end());
    xs.erase(unique(xs.begin(), xs.end()), xs.end());
    const int n_xs = xs.size();
    std::vector< kinStruct > ans(n_xs) ;
    // Get kernels in spec. For F2 there is only one value of t, i.e. 0
    std::vector<Reggeon> reggeons = spec[0].getReggeons();
    const int n_reggeons = reggeons.size();
    std::vector<double> kinematics, iznbars(n_reggeons,0);
    for(int i = 0; i < n_xs; i++)
    {
        kinematics = {xs[i]};
        for(int reg_index = 0; reg_index < n_reggeons; reg_index++) iznbars[reg_index] = IzNBar(kinematics, reggeons[reg_index], gs);
        kinStruct izbars(xs[i], iznbars);
        ans[i] = izbars;
    }
    return ans ;
}

std::vector<double> DeepInelasticScattering::predictObs(const std::vector<kinStruct> &Izs,
                                                        const std::vector<kinStruct> &IzsBar,
                                                        const std::vector< std::vector<double> > &points,
                                                        const bool savePredictions)
{
    // Check that Izs, IzsBar and points have the same size
    if (points.size() == 0) throw std::runtime_error("points argument has 0 length. Aborting DIS predictObs.");
    std::vector<double> ans(points[0].size()) ;
    kinStruct iznStruct, iznbarStruct;
    std::vector<double> izn, iznbar;
    for(int i = 0; i < points[0].size(); i++)
    {
        iznStruct = binary_search<kinStruct>(Izs, kinStruct(points[0][i],{}));
        iznbarStruct = binary_search<kinStruct>(IzsBar, kinStruct(points[1][i],{}));
        izn = iznStruct.izns; iznbar = iznbarStruct.izns;
        ans[i] = sum(izn * iznbar);
    }
    if(savePredictions)
    {
        std::ofstream myfile;
        std::string file_path;
        std::cout << "Please introduce the path to save the predictions of this DIS observable" << std::endl;
        std::cin >> file_path;
        myfile.open(file_path);
        myfile << "Q2\tx\tPred" << std::endl;
        for(int i = 0; i < points[0].size(); i++)
        {
            myfile << points[0][i] << '\t' << points[1][i] << '\t' << ans[i] << std::endl;
        }
    }
    return ans ;
}

DeepInelasticScattering& DeepInelasticScattering::operator=(const DeepInelasticScattering &rhs)
{
    if (this == &rhs) return *this;
    copy(rhs);
    return *this;
}

DeepInelasticScattering::~DeepInelasticScattering()
{}
                         

IzNIntegrand::IzNIntegrand(const Poly_Interp<double> &f1, const Poly_Interp<double> &f2, const U1NNMode &f3, const Poly_Interp<double> &f4):
                func1(f1), func2(f2), func3(f3), func4(f4) {}

double IzNIntegrand::operator()(const double x) {return 0;}

double f(double * x, void * params)
{
    return ((IzNIntegrand *) params)->operator()(*x);
}