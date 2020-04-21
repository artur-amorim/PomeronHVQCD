#ifndef ROOT_FIND_H
#define ROOT_FIND_H

#include <vector>
#include "../schrodinger/common.h"

template<typename T>
inline T SIGN(const T& a, const T &b)
{
	return (b) >= 0.0 ? std::fabs(a) : -std::fabs(a) ;
}

template<class T>
double zbrent(T& func, double x1, double x2, double tol, bool silent = false)
{
	const int ITMAX = 600;
	const double EPS = 1e-9;
	int iter;
	double a=x1,b=x2,c=x2,d,e,min1,min2;
	double fa=func(a),fb=func(b),fc,p,q,r,s,tol1,xm;

	if (!silent && ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)))
		std::cout << "?";
	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c  = a;			//Rename a, b, c and adjust bounding interval d.
			fc = fa;
			e  = d = b-a;
		}
		if (std::fabs(fc) < std::fabs(fb)) {
			a  = b;
			b  = c;
			c  = a;
			fa = fb;
			fb = fc;
			fc = fa;
		}

		tol1 = 2.0*EPS*std::fabs(b)+0.5*tol; //Convergence check.
		xm = 0.5*(c-b);

		if (std::fabs(xm) <= tol1 || fb == 0.0) return b;
		if (std::fabs(e) >= tol1 && std::fabs(fa) > std::fabs(fb)) {
			s=fb/fa;		//Attempt inverse quadratic interpolation.
			if (a == c) {
				p = 2.0*xm*s;
				q = 1.0-s;
			} else {
				q = fa/fc;
				r = fb/fc;
				p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q = (q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;		//Check whether in bounds.
			p=std::fabs(p);
			min1=3.0*xm*q-std::fabs(tol1*q);
			min2=std::fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e = d;
				//Accept interpolation.
				d = p/q;
			} else {
				d = xm;
				//Interpolation failed, use bisection.
				e = d;
			}
		} else {		//Bounds decreasing too slowly, use bisection.
			d = xm;
			e = d;
		}
		a = b;			//Move last best guess to a.
		fa = fb;
		if (std::fabs(d) > tol1)	//Evaluate new trial root.
			b += d;
		else
			b += SIGN(tol1,xm);
		fb = func(b);
	}
	std::cout << "[WARN] Maximum number of iterations exceeded in zbrent, returning biggest value"<< std::endl;
	return x2;
}

template<class T>
std::vector<Range> bracketZeros(T & func, const int n_zeros, const double x_start = 0, const double delta = 0.1)
{
    /*
        This function searches for the first n_zeros of the function func in steps of delta
        The first guessed interval is [0, delta] and is updated in steps of delta
        until we have n_zeros such ranges.
    */
    std::vector<Range> intervals;
    double x0 = x_start, x1 = x_start + delta;
	double f0 = func(x0), f1 = func(x1);
    int n = 0;
    // Start to search for zeros
    while(n < n_zeros)
    {
        if( f0 * f1 < 0)        // If there is a zero the function values at both ends have opposite signs
        {
            intervals.push_back(Range(x0, x1));     // append the interval
            x0 = x1;                                // Update x0 and x1
            x1 += delta;
			f0 = f1;
			f1 = func(x1);
            n++;
        }
        else
        {
            x0 = x1;                        // Didn't find zero, just update x0 and x1
            x1 += delta;
			f0 = f1;
			f1 = func(x1);
        }
    }
    return intervals;
}

#endif