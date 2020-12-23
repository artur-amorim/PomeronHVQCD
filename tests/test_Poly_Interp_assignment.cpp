#include <fstream>
#include <cmath>
#include "methods/interpolation/Poly_Interp.hpp"

using namespace std;


int main()
{
    vector<double> X, Y;
    for(double x = - 5; x <= 5; x+= 0.01)
    {
        X.push_back(x); Y.push_back(exp(-x*x));
    }
    Poly_Interp<double> f(X, Y, 4);
    Poly_Interp<double> f_copy;
    f_copy = f;

    for(double x = -3.2123; x <= 4.12313; x += 0.12314) cout << x << '\t' << f_copy.interp(x) << '\t' << exp(-x*x) << endl;

    return 0;
}