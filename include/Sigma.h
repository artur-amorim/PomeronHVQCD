#ifndef SIGMA_H
#define SIGMA_H

#include <vector>
#include <string>
#include "Process.h"
#include "Reggeon.h"

class Sigma: public Process
{
    private:
        double barn_to_GEVMINUS2;
        void loadData(const std::string & file_path, const double conv_factor);
        void copy(const Sigma &sigma);
    public:
        Sigma(const std::string &file_path, const double conv_factor);
        Sigma(const Sigma &sigma);
        std::vector<double> expVal();
        std::vector<double>  expErr();
        std::vector<std::vector<double> >  expKinematics();
        virtual double IzN(const std::vector<double> &kin, const Reggeon &reg);
        double IzNBar(const std::vector<double> &kin, const Reggeon &reg, const std::vector<double> &gs);
        std::vector<kinStruct>  getIzs(const std::vector< std::vector<double> > &points, const std::vector<Spectra> &spec);
        std::vector<kinStruct>  getIzsBar(const std::vector< std::vector<double> >  &points, const std::vector<Spectra> &spec,
                                                  const std::vector<double>  &gs);
        std::vector<double>  predict(const std::vector<kinStruct>  &Izs, const std::vector<kinStruct>  &IzsBar,
                                            const std::vector< std::vector<double> >  &points, const bool savePredictions = false);
        std::vector<double>  diffObsWeighted(const std::vector<kinStruct>  &Izs, const std::vector<kinStruct>  &IzsBar, 
                                                    const std::vector< std::vector<double> >  &points);
        ~Sigma();
};

#endif