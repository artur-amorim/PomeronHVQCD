#ifndef NELDERMEAD_HPP
#define NELDERMEAD_HPP

#include <vector>
#include <algorithm>
#include <exception>
#include <armadillo>
#include <cmath>

// The following code was adapted from Numerical Recipes and from Gao et. al
// Implementing the Nelder-Mead simplex algorithm with adaptive parameters
struct NelderMead {
    // Minimization of a muldimensional function by the Nelder-Mead method
	const double ftol;
	int nfuncEval;                          // Number of function evaluations
	int mpts;
	int ndim;
	double minf;                        // Minimum value of the function
	std::vector<double> y;                      // Function values at the simplex
	arma::mat p;                      // Simplex
	// Constructor - ftoll is the fractional convergence tolerance to be achieved in the function value
    NelderMead(const double ftoll) : ftol(ftoll) {}

	template <class T>
	std::vector<double> minimize(std::vector<double> &point, const double del, T &func)
	{
		/*
			Multidimensional minimization of the function or functor func(x) by the simplex method
			of Nelder and Mead. x = {x0, x1, ... x_{ndim - 1}}.
			Returns the position of the minimum
		*/
		std::vector<double> dels(point.size(),del);
		return minimize(point,dels,func);
	}

	template <class T>
	std::vector<double> minimize(std::vector<double> &point, const std::vector<double> &dels, T &func)
	{
        /*
			Alternative interface that takes different displacements 
			dels = {del0, del1, ... del(ndim-1)} in different directions
			for the initial simplex.
        */
		int ndim=point.size();
		arma::mat simplex(ndim+1,ndim);
		for (int i=0;i<ndim+1;i++) {
			for (int j=0;j<ndim;j++)
				simplex(i,j)=point[j];
			if (i !=0 ) simplex(i,i-1) += dels[i-1];
		}
		return minimize(simplex,func);
	}

	template <class T>
	std::vector<double> minimize(arma::mat &simplex, T &func)
	{
		/*
			Most general interface: initial simplex specified by the matrix simplex[0..ndim][0..ndim-1].
        	Its ndim+1 rows are ndim-dimensional vectors that are the vertices of the starting simplex.
		*/
        ndim = simplex.n_cols;
		const int maxfuncEval= 1000 * ndim;		// Maximum allowed number of function evaluations.
		const double epsilon = 1.0e-9;
		int ihi,ilo,inhi;
		mpts=simplex.n_rows;
		ndim=simplex.n_cols;
		const double alpha = 1;
		const double beta = 1 + 2.0 / ndim;
		const double gamma_fac = 0.75 - 0.5 / ndim ;
		const double delta = 1 - 1.0 / ndim;
		std::vector<double> psum(ndim),pmin(ndim),x(ndim);
		p=simplex;
		y.resize(mpts);
		for (int i=0;i<mpts;i++) {
			for (int j=0;j<ndim;j++)
				x[j]=p(i,j);
			y[i]=func(x);
		}
		nfuncEval=0;
		get_psum(p,psum);
		for (;;) {
			ilo=0;
			/*
				First we must determine which point is the highest(worst), next-highest
				and lowest (best), by looping over the points in the simples.
			*/
			ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
			for (int i=0;i<mpts;i++) {
				if (y[i] <= y[ilo]) ilo=i;
				if (y[i] > y[ihi]) {
					inhi=ihi;
					ihi=i;
				} else if (y[i] > y[inhi] && i != ihi) inhi=i;
			}
			double rtol=2.0*std::fabs(y[ihi]-y[ilo])/(std::fabs(y[ihi])+std::fabs(y[ilo])+epsilon);
			// Coompute the fractional range from highest to lowest and return if satisfactory.
			if (rtol < ftol) {
				std::swap(y[0],y[ilo]);			// If returning, put best point and value in slot 0
				for (int i=0;i<ndim;i++) {
					std::swap(p(0,i),p(ilo,i));
					pmin[i]=p(0,i);
				}
				minf=y[0];
				return pmin;
			}
			if (nfuncEval >= maxfuncEval) throw std::runtime_error("Maximum function evaluations exceeded");
			nfuncEval += 2;
			/*
				Begin a new iteration. First extrapolate by a factor -1 through the face of the simplex
				across from the high point, i.e. reflect the simplex from the high point.
			*/
			double ytry=simptry(p,y,psum,ihi,-alpha,func);
			if (ytry <= y[ilo])
				ytry=simptry(p,y,psum,ihi,beta,func);	// Gives a result better than the best point, so try an additional extrapolation by factor 2.
			else if (ytry >= y[inhi]) {
				double ysave=y[ihi];				// The reflected point is worse than the second-highest, so look for an intermediate lower point, i.e. one-dim contraction.
				ytry=simptry(p,y,psum,ihi,gamma_fac,func);
				if (ytry >= ysave) {		// Can't seem to get rid of that high point.
					for (int i=0;i<mpts;i++) {		// Better contract around the lowest (best) point.
						if (i != ilo) {
							for (int j=0;j<ndim;j++)
								p(i,j)=psum[j]= delta * p(i,j) + (1-delta) * p(ilo,j);
							y[i]=func(psum);
						}
					}
					nfuncEval += ndim;			// Keep track of function evaluations
					get_psum(p,psum);			// Recompute psum
				}
			} else --nfuncEval;					// Correct the evaluation count.
		}										// Go back for the test of doneness and the next iteration.
	}
	inline void get_psum(arma::mat &p, std::vector<double> &psum)
	{
		// Utility function
		for (int j=0;j<ndim;j++) {
			double sum=0.0;
			for (int i=0;i<mpts;i++)
				sum += p(i,j);
			psum[j]=sum;
		}
	}

	template <class T>
	double simptry(arma::mat &p, std::vector<double> &y, std::vector<double> &psum,
		const int ihi, const double fac, T &func)
	{
		/*
			Helper function: Extrapolates by a factor fac through the face of the simplex
			across from the high point, tries it, and replaces the high point
			if the new point is better.
		*/
		std::vector<double> ptry(ndim);
		double fac1=(1.0-fac)/ndim;
		double fac2=fac1-fac;
		for (int j=0;j<ndim;j++)
			ptry[j]=psum[j]*fac1-p(ihi,j)*fac2;
		double ytry=func(ptry);			// Evaluate the function at the trial point

		if (ytry < y[ihi]) {			// If it's better than the highest, then replace the highest.
			y[ihi]=ytry;
			for (int j=0;j<ndim;j++) {
				psum[j] += ptry[j]-p(ihi,j);
				p(ihi,j)=ptry[j];
			}
		}
		return ytry;
	}
};


template<class T>
std::vector<double> optimFunction(const std::vector<double> &x, T &func, const double delta, const double ftol = 1e-3)
{
    // Given an initial guess of parameters X,
    // prints the value for which the function or functor f has a minimum
    NelderMead NM(ftol) ;
	// Start the optimization process
	std::vector<double> xguess = x;
	std::vector<double> optParams = NM.minimize(xguess, delta, func);
    return optParams;
};

template<class T>
std::vector<double> optimFunction(const std::vector<double> &x, T &func, const std::vector<double> &deltas, const double ftol = 1e-3) 
{
    // Given an initial guess of parameters X,
    // prints the value for which the function or functor f has a minimum
    NelderMead NM(ftol) ;
	// Start the optimization process
	std::vector<double> xguess = x;
	std::vector<double> optParams = NM.minimize(xguess, deltas, func);
	// Return the optimal set of parameters
    return optParams;
};

#endif