#ifndef NELDERMEAD_HPP
#define NELDERMEAD_HPP

#include <vector>
#include <algorithm>
#include <exception>
#include <armadillo>
#include <cmath>

// The following code was adapted from Numerical Recipes
struct NelderMead {
    // Minimization of a muldimensional function by the Nelder-Mead method
	const double ftol;
	int nfuncEval;                          // Number of function evaluations
	int mpts;
	int ndim;
	double minf;                        // Minimum value of the function
	std::vector<double> y;                      // Function values at the simplex
	arma::mat p;                      // Simplex
	// Constructor 
    NelderMead(const double ftoll) : ftol(ftoll) {}

	template <class T>
	std::vector<double> minimize(std::vector<double> &point, const double del, T &func)
	{
		std::vector<double> dels(point.size(),del);
		return minimize(point,dels,func);
	}

	template <class T>
	std::vector<double> minimize(std::vector<double> &point, std::vector<double> &dels, T &func)
	{
        // Minimization of the function or functor func(x), where x has dimension n
        // The method used is the Nelder and Mead one
        // The initial simplex is provided in point and the displacement along each direction in dels
        // This routine returns the location of the minimum
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
		// Most general interface: initial simplex specified by the matrix simplex[0..ndim][0..ndim-1].
        // Its ndim+1 rows are ndim-dimensional vectors that are the vertices of the starting simplex.
        const int maxfuncEval=10000;
		const double epsilon=1.0e-9;
		int ihi,ilo,inhi;
		mpts=simplex.n_rows;
		ndim=simplex.n_cols;
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
			ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
			for (int i=0;i<mpts;i++) {
				if (y[i] <= y[ilo]) ilo=i;
				if (y[i] > y[ihi]) {
					inhi=ihi;
					ihi=i;
				} else if (y[i] > y[inhi] && i != ihi) inhi=i;
			}
			double rtol=2.0*std::fabs(y[ihi]-y[ilo])/(std::fabs(y[ihi])+std::fabs(y[ilo])+epsilon);
			if (rtol < ftol) {
				std::swap(y[0],y[ilo]);
				for (int i=0;i<ndim;i++) {
					std::swap(p(0,i),p(ilo,i));
					pmin[i]=p(0,i);
				}
				minf=y[0];
				return pmin;
			}
			if (nfuncEval >= maxfuncEval) throw std::runtime_error("Maximum function evaluations exceeded");
			nfuncEval += 2;
			double ytry=simptry(p,y,psum,ihi,-1.0,func);
			if (ytry <= y[ilo])
				ytry=simptry(p,y,psum,ihi,2.0,func);
			else if (ytry >= y[inhi]) {
				double ysave=y[ihi];
				ytry=simptry(p,y,psum,ihi,0.5,func);
				if (ytry >= ysave) {
					for (int i=0;i<mpts;i++) {
						if (i != ilo) {
							for (int j=0;j<ndim;j++)
								p(i,j)=psum[j]=0.5*(p(i,j)+p(ilo,j));
							y[i]=func(psum);
						}
					}
					nfuncEval += ndim;
					get_psum(p,psum);
				}
			} else --nfuncEval;
		}
	}
	inline void get_psum(arma::mat &p, std::vector<double> &psum)
	{
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
		std::vector<double> ptry(ndim);
		double fac1=(1.0-fac)/ndim;
		double fac2=fac1-fac;
		for (int j=0;j<ndim;j++)
			ptry[j]=psum[j]*fac1-p(ihi,j)*fac2;
		double ytry=func(ptry);
		if (ytry < y[ihi]) {
			y[ihi]=ytry;
			for (int j=0;j<ndim;j++) {
				psum[j] += ptry[j]-p(ihi,j);
				p(ihi,j)=ptry[j];
			}
		}
		return ytry;
	}
};

struct Fit
{
	std::vector<double> bestFitPars;
	std::vector<double> bestFitParsUncert;
	Fit(const std::vector<double> &pars, const std::vector<double> &uncert):
		bestFitPars(pars), bestFitParsUncert(uncert) {}
};

std::vector<double> Pij(std::vector<double> &Pi, std::vector<double> &Pj)
{
	if (Pi.size() != Pj.size()) throw std::runtime_error("Vectors Pi and Pj with different sizes in Pij function in NelderMead.hpp");
	const int n = Pi.size();
	std::vector<double> ans(n,0);
	// Perform the sum
	for(int k = 0; k < n; k++) ans[k] = 0.5 * (Pi[k] + Pj[k]);
	return ans;
}

template<class T>
std::vector<double> parameterUncertainties(T &func, const NelderMead &NM)
{
	// Computes the parameter uncertainties
	// Get the simplex and its function values
	arma::mat simp = NM.p;
	const int nSimPts = simp.n_rows;
	const int n = simp.n_cols;
	std::vector<std::vector<double> > simplex;
	for(int i = 0; i < nSimPts; i++)
	{
		std::vector<double> row(n,0);
		for (int j = 0; j < n; j++) row[j] = simp(i,j);
		simplex.push_back(row);
	}
	std::vector<double> simplexFuncVals = NM.y;
	// Compute B matrix
	arma::mat B(n, n);
	for(int i = 1; i <= n; i++)
	{
		for(int j = 1; j <= n; j++)
		{
			if (i > j) B(i-1,j-1) = B(j-1,i-1);
			else if (i == j)
			{
				double yi = simplexFuncVals[i];
				double y0 = simplexFuncVals[0];
				std::vector<double> p0i = Pij(simplex[0], simplex[i]);
				double yi0 = func(p0i);
				B(i-1,i-1) = 2 * (yi + y0 - 2 * yi0);
			}
			else
			{
				double y0 = simplexFuncVals[0];
				std::vector<double> pij = Pij(simplex[i], simplex[j]);
				double yij = func(pij);
				std::vector<double> p0i = Pij(simplex[0], simplex[i]);
				double y0i = func(p0i);
				std::vector<double> p0j = Pij(simplex[0], simplex[j]);
				double y0j = func(p0j);
				B(i-1, j-1) = 2 * (yij + y0 - y0i - y0j);
			}
		}
	}
	// Compute Q matrix
	arma::mat Q(n, n);
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++) Q(i,j) = simplex[j+1][i] - simplex[0][i];
	}
	// We can now build the covariance matrix from B and Q
	arma::mat CovMat = Q * arma::inv(B) * Q.t();
	// The parameters uncertainties are the square roots of the diagonal elements of CovMat
	std::vector<double> parUncert(n, 0);
	for(int i = 0; i < n; i++) parUncert[i] = std::sqrt(CovMat(i,i));
	return parUncert;
}

template<class T>
Fit optimFunction(const std::vector<double> &x, T &func, double delta)
{
    // Given an initial guess of parameters X,
    // prints the value for which the function or functor f has a minimum
    NelderMead NM(1e-3) ;
	// Start the optimization process
	std::vector<double> xguess = x;
	std::vector<double> optParams = NM.minimize(xguess, delta, func);
	// Compute the parameterUncertainties
	std::vector<double> parUncert = parameterUncertainties(func, NM);
    return Fit(optParams, parUncert);
};

template<class T>
Fit optimFunction(const std::vector<double> &x, T &func, std::vector<double> deltas)
{
    // Given an initial guess of parameters X,
    // prints the value for which the function or functor f has a minimum
    NelderMead NM(1e-3) ;
	// Start the optimization process
	std::vector<double> xguess = x;
	std::vector<double> optParams = NM.minimize(xguess, deltas, func);
	// Compute the parameterUncertainties
	std::vector<double> parUncert = parameterUncertainties(func, NM);
	// Return the optimal set of parameters
    return Fit(optParams, parUncert);
};

#endif