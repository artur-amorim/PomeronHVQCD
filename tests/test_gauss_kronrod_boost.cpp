#include <iostream>
#include <iomanip>
#include <cmath>
#include <boost/math/quadrature/gauss_kronrod.hpp>

using namespace std;
using namespace boost::math::quadrature;

class func
{
    private:
        double alpha;
    public:
        func(const double a): alpha(a) {}
        double operator() (const double x) const { return exp(-x*x / alpha) ;}
};

int main()
{
    func f(3);
    double error;
    double Q = gauss_kronrod<double, 61>::integrate(f, 0, std::numeric_limits<double>::infinity(), 15, 1e-15, &error);

    cout << fixed;
    cout << setprecision(15);

    cout << "Q = " << Q << " +-" << error << endl;

    return 0;
}