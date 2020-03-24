#include <iostream>
#include <vector>
#include "methods/optimization/NelderMead.hpp"

using namespace std;

static int f_calls(0);

double f(vector<double> &x)
{
    f_calls += 1;
    double sum = 0;
    for(int i = 0; i < x.size(); i++) sum += x[i] * x[i]; 
    return sum;
}

int main(int argc, char ** argv)
{
    int ndim = 0;
    if (argc < 2) ndim = 2;
    else ndim = stoi(argv[1]);
    vector<double> x_guess(ndim, 1);
    vector<double> x_opt = optimFunction(x_guess, f, 0.1, 1e-3);
    cout << "Number of function evaluations: " << f_calls << endl;
    for(int i = 0; i < x_opt.size(); i++) cout << x_opt[i] << '\t';
    cout << endl;
    return 0;
}