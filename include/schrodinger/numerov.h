#ifndef NUMEROV_H
#define NUMEROV_H

#include "common.h"
#include "solvspec.h"

class Numerov : public SolvSpec{
    private:
        std::vector<Range> zeros;
        std::vector<Point> solLR, solRL, sol;
        Point minPot;
        void findMinPot();
        void buildSol();
        double diff(double E);
        void scanForZeroRanges(int nZeros);
        int getSolIndex();
        std::function<double(double)> diffFunc = [this](double E){ return diff(E);};
    public:
        Numerov();
        void findSpectrum(int nEigen);
        ~Numerov();
};
;
#endif