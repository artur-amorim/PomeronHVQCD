#ifndef YANGMILLS_H
#define YANGMILLS_H

#include "background.h"

class YangMills: public Background
{
    private:
        // Yang Mills EOMs
        void eom(const state &X, state &dXdA, const double A);
        // defines how to save the computed background fields for Yang-Mills
        void observer(const state &X, const double A);
        // prepare background for further computations
        void finalizeBackground();
    public:
        // YangMills constructor
        YangMills(const double ssc = 3.0, const double VVgIR = 2.05);
        // YangMills copy constructor
        YangMills(const YangMills &ym);
        // YangMills assignment operator
        YangMills& operator=(const YangMills &rhs); 
        // YangMills destructor
        ~YangMills();
        // Solves Yang Mills equations of motion
        void solve();
};

#endif