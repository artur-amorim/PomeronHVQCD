#include <iostream>
#include "HolographicVQCD.h"
#include "U1NNMode.h"

using namespace std;

int main(int argc, char ** argv)
{
    double Q2 = stod(argv[1]);

    U1NNMode mode(Q2);
    mode.computeMode();
    cout << "Mode computed" << endl;
    mode.saveMode("plots/U1NNMode/U1NNMode_Q2_"+to_string(Q2)+".txt");
    return 0;
}