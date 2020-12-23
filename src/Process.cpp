#include "HolographicVQCD.h"
#include "Process.h"
#include "methods/vectorOperators.hpp"

// Constructor of kinStruct
kinStruct::kinStruct(const double var, const std::vector<double> &izs): kinVar(var), izns(izs) {}

bool kinStruct::operator<(const kinStruct &rhs) const
{
    return kinVar < rhs.kinVar ;
}

bool kinStruct::operator>(const kinStruct &rhs) const
{
    return kinVar > rhs.kinVar ;
}

std::vector<double> Process::u = {};
std::vector<double> Process::z = {};
std::vector<double> Process::Astring = {};
std::vector<double> Process::GluonPotFac = {};
std::vector<double> Process::MesonPotFac = {};
Poly_Interp<double> Process::GluonPotFactor({},{},4);
Poly_Interp<double> Process::MesonPotFactor({}, {}, 4);
Poly_Interp<double> Process::ufunc({},{},4);


void Process::copy(const Process &rhs)
{
    dataPoints = rhs.dataPoints;
}

Process::Process(): dataPoints({})
{
    // We now compute potFactor
    if (z.size() == 0)
    {
        u = hvqcd().u(); z = hvqcd().z();
        Astring = hvqcd().Astring();
        std::vector<double> Gs = hvqcd().G(); 
        GluonPotFac.resize(z.size()); MesonPotFac.resize(z.size());
        // e^(-10/3 \Phi) V_f w_s^2 = e^(-2 \Phi / 3) V_f w^2 in terms of the Einstein frame background potentials
        std::vector<double> Phis = hvqcd().Phi(), taus = hvqcd().tau();
        for(int i = 0; i < z.size(); i++) 
        {
            GluonPotFac[i] = std::exp(-0.5 * Astring[i]);
            MesonPotFac[i] = std::sqrt(std::exp(-2*Phis[i]/3) * hvqcd().Vf(Phis[i], taus[i]) * std::pow(hvqcd().w(Phis[i]),2));
        }
        // Now we reverse because we want them from the UV to the IR
        std::reverse(std::begin(u), std::end(u));
        std::reverse(std::begin(z), std::end(z));
        std::reverse(std::begin(Astring), std::end(Astring));
        std::reverse(std::begin(GluonPotFac), std::end(GluonPotFac));
        std::reverse(std::begin(MesonPotFac), std::end(MesonPotFac));
        ufunc = Poly_Interp<double>(z, u, 4);
        GluonPotFactor = Poly_Interp<double>(z, GluonPotFac, 4);
        MesonPotFactor = Poly_Interp<double>(u, MesonPotFac, 4);
    }
}

Process::Process(const Process &proc): dataPoints(proc.dataPoints){}

std::vector<std::vector<double> > Process::getDataPts()
{
    return dataPoints;
}

void Process::setDataPts(const std::vector<std::vector<double> > &pts)
{
    dataPoints = pts ;
}

std::vector<double>  Process::getNeededTVals()
{
    std::vector<double> pts = {0};
    return pts ;
}

double Process::chi2(const std::vector<kinStruct> &Izs, const std::vector<kinStruct> &IzsBar, const std::vector< std::vector < double > > &points)
{
    const std::vector<double> obs = this->diffObsWeighted(Izs, IzsBar, points);
    const int n = obs.size() ;
    double ans = 0.0 ;
    for(int i = 0; i < n ; i++) ans += obs[i] * obs[i] ;
    return ans ;
}

Process& Process::operator= (const Process &proc)
{
    if (this == &proc) return *this;
    copy(proc);
    return *this;
}

Process::~Process(){}