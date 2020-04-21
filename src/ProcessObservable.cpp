#include "ProcessObservable.h"
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

ProcessObservable::ProcessObservable(const bool rrsslog): rsslog(rrsslog), datapts({}){}

ProcessObservable::ProcessObservable(const ProcessObservable &proc): rsslog(proc.rsslog), datapts(proc.datapts){}

bool ProcessObservable::getRSSLOG()
{
    return rsslog;
}

std::vector<std::vector<double> > ProcessObservable::getDataPts()
{
    return datapts;
}

void ProcessObservable::setRSSLOG(const bool rrsslog)
{
    rsslog = rrsslog ;
}

void ProcessObservable::setDataPts(const std::vector<std::vector<double> > &pts)
{
    datapts = pts ;
}

void ProcessObservable::copy(const ProcessObservable &rhs)
{
    rsslog = rhs.rsslog;
    datapts = rhs.datapts;
}

void  ProcessObservable::loadData(std::string file_path){}


std::vector<double> ProcessObservable::expVal() 
{
    std::vector<double> pts(0);
    return pts ;
}


std::vector<double> ProcessObservable::expErr()
{
    std::vector<double> pts(0);
    return pts ;
}


std::vector< std::vector<double> >  ProcessObservable::expKinematics()
{
    std::vector< std::vector<double> > pts(0) ;
    return pts ;
}


std::vector<double>  ProcessObservable::getNeededTVals()
{
    std::vector<double> pts(0);
    return pts ;
}

double ProcessObservable::IzN(const std::vector<double> &kin, const Reggeon &reg)
{
    return 0.0 ;
}

double ProcessObservable::IzNBar(const std::vector<double> &kin, const Reggeon &reg, const std::vector<double> &gs)
{
    return 0.0 ;
}


std::vector<kinStruct> ProcessObservable::getIzs(const std::vector< std::vector<double> >  &points, const std::vector<Spectra> &spec)
{
    std::vector<kinStruct> pts;
    return pts ;
}


std::vector<kinStruct>  ProcessObservable::getIzsBar(const std::vector< std::vector<double> >  &points,
                                                     const std::vector<Spectra> &spec,
                                                     const std::vector<double>  &gs)
{
    std::vector<kinStruct> pts;
    return pts ;
}


std::vector<double>  ProcessObservable::predictObs(const std::vector<kinStruct>  &Izs, const std::vector<kinStruct>  &IzsBar,
                                                  const std::vector< std::vector<double> >  &points, const bool savePredictions)
{
    std::vector<double> pts(0) ;
    return pts ;
}


std::vector<double> ProcessObservable::diffObsWeighted(const std::vector<kinStruct> &Izs, const std::vector<kinStruct> &IzsBar, const std::vector< std::vector < double > > &points)
{
    if( points.size() == 0) std::vector< std::vector< double > > points = this->expKinematics() ;   // If points is NULL provide the experimental ones
    const std::vector<double> Opred = this->predictObs(Izs, IzsBar, points, false) ;                // Predictions of the model for the process
    const std::vector<double> Oexp  = this->expVal() ;                                              // Experimental values of the process
    const std::vector<double> Oerr  = this->expErr() ;                                              // Experimental errors of the process
    if (rsslog) return Oexp * log(Opred / Oexp) / Oerr ;
    else return (Opred - Oexp) / Oerr ;

}

double ProcessObservable::rss(const std::vector<kinStruct> &Izs, const std::vector<kinStruct> &IzsBar, const std::vector< std::vector < double > > &points)
{
    const std::vector<double> obs = this->diffObsWeighted(Izs, IzsBar, points);
    const int n = obs.size() ;
    double chi2 = 0.0 ;
    for(int i = 0; i < n ; i++) chi2 += obs[i] * obs[i] ;
    return chi2 ;
}

ProcessObservable& ProcessObservable::operator= (const ProcessObservable &proc)
{
    if (this == &proc) return *this;
    copy(proc);
    return *this;
}

ProcessObservable::~ProcessObservable(){}