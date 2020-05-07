#ifndef MATRIX_NUMEROV_H
#define MATRIX_NUMEROV_H

#include "solvspec.h"

class MatrixNumerov : public SolvSpec{
    public:
        MatrixNumerov();
        void findSpectrum(const int nEigen);
        ~MatrixNumerov();
};
;
#endif