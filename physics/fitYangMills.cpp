#include <iostream>
#include <vector>
#include "YangMills.h"
#include "schrodinger/chebspec.h"
#include "schrodinger/schrodinger.h"

using namespace std;

int main(int argc, char ** argv)
{
    double sc_guess = 3, VgIR_guess = 2.05;
    if (argc == 3)
    {
        sc_guess = std::stod(argv[1]); VgIR_guess = std::stod(argv[2]);
    }

    // Set number of Chebyshev points 
    chebSetN(800);

    // Solve Yang-Mills for sc = 3, VgIR = 2.05
    YangMills ym(sc_guess, VgIR_guess);
    ym.solve(-80, 20);

    // Compute the initial mass ratios
    computeYangMillsRatios(ym);
    
    std::cout << "Now fitting" << endl;
    fitYangMills(sc_guess, VgIR_guess);

    return 0;
}