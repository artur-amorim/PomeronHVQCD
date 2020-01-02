#ifndef POLY_INTERP_HPP
#define POLY_INTERP_HPP

#include <iostream>
#include "Base_Interp.hpp"

template<class T>
class Poly_Interp: public Base_Interp<T>
{
	// Polynomial interpolation object. Construct with x and y vectors, and the number M of points
	// to be used locally (polynomial order plus one), then call interp for interpolated values
	private:
		double dy;			// Error estimate dy for the most recent call
		double rawinterp(int jlo, T x);
		void copy(const Poly_Interp<T> &poly);
	public:
		Poly_Interp();
		Poly_Interp(const std::vector<T> &X, const std::vector<T> &Y, int m);
		Poly_Interp(const Poly_Interp<T> &poly);
		double getDY();		// Returns dy
		Poly_Interp<T>& operator= (const Poly_Interp<T> &poly);
		~Poly_Interp();
};

template<class T>
Poly_Interp<T>::Poly_Interp(): Base_Interp<T>(), dy(0.0) {}

template<class T>
Poly_Interp<T>::Poly_Interp(const std::vector<T> &X, const std::vector<T> &Y, int m): Base_Interp<T>(X, Y, m), dy(0.0) {}

template<class T>
void Poly_Interp<T>::copy(const Poly_Interp<T> &poly)
{
	Base_Interp<T>::copy(poly);
	dy = poly.dy;
}

template<class T>
Poly_Interp<T>::Poly_Interp(const Poly_Interp<T> &poly):
	Base_Interp<T>::Base_Interp(poly), dy(poly.dy) {}

template<class T>
Poly_Interp<T>& Poly_Interp<T>::operator= (const Poly_Interp<T> &poly)
{
	// We don't want to waste resources 
	if (this == &poly) {return *this;}
	copy(poly);
	return *this;
}

template<class T>
double Poly_Interp<T>::getDY()
{
	return this->dy;
}

template<class T>
double Poly_Interp<T>::rawinterp(int jl, T x)
{
	/*
		Given a value x, and using pointers to data X and Y, this routine returns an interpolated value y, and stores an error estimate dy.
		The returned value is obtained by mm-point polynomial interpolation on the subrange X[jl..jl+mm-1]
	*/
	int m = this->getM();
	int i, mm, ns = 0;
	double y, den, dif, dift, ho, hp, w;
	std::vector<double> c(m), d(m);
	dif = std::fabs(x-this->X(jl));
	for(i = 0; i < m; i++)
	{
		if((dift=std::fabs(x-this->X(i+jl))) < dif)
		{
			ns = i;
			dif = dift;
		}
		c[i] = this->Y(i+jl);
		d[i] = this->Y(i+jl);
	}
	y = this->Y(ns-- + jl);
	for(mm = 1; mm < m; mm++)
	{
		for(i = 0; i < m - mm; i++)
		{
			ho = this->X(i+jl) - x;
			hp = this->X(i+mm+jl) - x;
			w = c[i+1] - d[i];
			if ((den=ho-hp) == 0.0) throw std::runtime_error("Poly_Interp::rawinterp division by zero");
			den = w /den;
			d[i] = hp * den;
			c[i] = ho * den;
		}
		y += (dy=(2*(ns+1) < (m - mm) ? c[ns+1] : d[ns--]));
	}
	return y;
}

template<class T>
Poly_Interp<T>::~Poly_Interp() {}

#endif