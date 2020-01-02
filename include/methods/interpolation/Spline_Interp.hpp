#ifndef SPLINE_INTERP_HPP
#define SPLINE_INTERP_HPP

#include <iostream>
#include "Base_Interp.hpp"

template<class T>
class Spline_Interp:public Base_Interp<T>
{
	private:
		std::vector<T> Y2 ; // Values of the 2nd derivatives
		void sety2(const std::vector<T> &X, const std::vector<T> &Y, T yp1 = 1e99, T ypn = 1e99);
		T rawinterp(int jl, T x);
        T rawder1(int jl, T x);
        T rawder2(int jl, T x);
		void copy(const Spline_Interp<T> &spline);
	public:
		Spline_Interp();
		Spline_Interp(const std::vector<T> &X, const std::vector<T> &Y, T yp1 = 1e99, T ypn = 1e99);
		Spline_Interp(const Spline_Interp<T> &spline);
		std::vector<T> getY2() const;
		T integrate (T x);
		Spline_Interp<T>& operator= (const Spline_Interp<T> & spline);
		~Spline_Interp();
};

template<class T>
Spline_Interp<T>::Spline_Interp(): Base_Interp<T>(), Y2({}) {}

template<class T>
Spline_Interp<T>::Spline_Interp(const std::vector<T> &X, const std::vector<T> &Y, T yp1, T ypn) : Base_Interp<T>(X, Y, 2), Y2({})
{
	Y2.resize(X.size()) ;
	sety2(X, Y, yp1, ypn); // Computes the values of y2
}

template<class T>
void Spline_Interp<T>::sety2(const std::vector<T> &X, const std::vector<T> &Y, T yp1, T ypn)
{
	T p, qn, sig, un ;
	int n = Y2.size() ;
	std::vector<T> U(n-1) ;
	if (yp1 > 0.99e99)
	{
		Y2[0] = U[0] = 0.0 ;
	}
	else
	{
		Y2[0] = - 0.5 ;
		U[0] = (3.0/(X[1]-X[0])) * ((Y[1]-Y[0])/(X[1]-X[0])-yp1);
	}
	for(int i = 1; i < n - 1; i++)
	{
		sig = (X[i]-X[i-1]) / (X[i+1]-X[i-1]) ;
		p = sig * Y2[i-1] + 2.0 ;
		Y2[i] = (sig -1.0) / p ;
		U[i] = (Y[i+1]-Y[i]) / (X[i+1]-X[i]) - (Y[i]-Y[i-1]) / (X[i] - X[i-1]) ;
		U[i] = (6.0 * U[i] / (X[i+1]-X[i-1]) - sig * U[i-1]) / p ;
	}
	if (ypn > 0.99e99)
	{
		qn = un = 0.0 ;
	}
	else
	{
		qn = 0.5 ;
		un = (3.0 / (X[n-1] - X[n-2])) * (ypn - (Y[n-1] - Y[n-2]) / (X[n-1] - X[n-2])) ;
	}
	Y2[n-1] = (un - qn * U[n-2]) / (qn * Y2[n-2] + 1.0 );
	for (int k = n - 2; k >= 0 ; k--)
	{
		Y2[k] = Y2[k] * Y2[k+1] + U[k];
	}
}

template<class T>
Spline_Interp<T>::Spline_Interp(const Spline_Interp<T> &spline):
	Base_Interp<T>::Base_Interp(spline), Y2(spline.Y2) {}

template<class T>
T Spline_Interp<T>::rawinterp(int jl, T x)
{
	int klo = jl, khi = jl + 1 ;
	T y, h, b, a ;
	h = this->X(khi) - this->X(klo) ;
	if ( h == 0.0) throw std::runtime_error("Spline_Interp::rawinterp division by zero");
	a = ( this->X(khi) - x ) / h ;
	b = ( x - this->X(klo) ) / h ;
	y = a * this->Y(klo) + b * this->Y(khi) + ( (a*a*a - a) * Y2[klo] + (b*b*b - b) * Y2[khi] ) * ( h * h) / 6.0;
	return y ;
}

template<class T>
T Spline_Interp<T>::rawder1(int jl, T x)
{
	int klo = jl, khi = jl + 1 ;
	T dy, h, b, a ;
	h = this->X(khi) - this->X(klo) ;
	if ( h == 0.0) throw std::runtime_error("Spline_Interp::rawder1 division by zero");
	a = ( this->X(khi) - x ) / h ;
	b = ( x - this->X(klo) ) / h ;
	dy = (this->Y(khi) - this->Y(klo))/h - (3*a*a - 1) * h * Y2[klo]/6.0 + (3*b*b - 1) * h * Y2[khi] / 6.0 ;
	return dy ;
}

template<class T>
T Spline_Interp<T>::rawder2(int jl, T x)
{
	int klo = jl, khi = jl + 1 ;
	T d2y, h, b, a ;
	h = this->X(khi) - this->X(klo) ;
	if ( h == 0.0) throw std::runtime_error("Spline_Interp::rawder2 division by zero");
	a = ( this->X(khi) - x ) / h ;
	b = ( x - this->X(klo) ) / h ;
	d2y = a * Y2[klo] + b * Y2[khi] ;
	return d2y ;
}

template<class T>
T Spline_Interp<T>::integrate (T x)
{
	if (x < this->X(0) || x > this->X(this->getN() - 1)) throw std::runtime_error("Spline_Interp::integrate Value out of bounds");
	int j = this->getCor() ? this->hunt(x) : this->locate(x) ; // Find j for x_j <= x < x_(j+1)
	T value = 0 ;
	// See notebook for formula of int_{x_j}^{x_(j+1)} f(x) dx where f is a cubic spline
	for (int k = 0 ; k < j ; k++)
	{
		T xk = this->X(k) ;
		T xk1 = this->X(k+1) ;
		value += ( xk - xk1 ) * ( -12*this->Y(k) -12*this->Y(k+ 1) + std::pow(xk - xk1, 2.0) * ( Y2[j] + Y2[j + 1] ) ) / 24.0 ;
	}
	T xj = this->X(j) ;
	T xj1 = this->X(j+1) ;
	value += ( x - xj ) * ( 12 * ( x + xj - 2 * xj1 ) * this->Y(j) +
	( x - xj ) * ( -12*this->Y(j+1) + std::pow(x + xj - 2 * xj1, 2.0) * Y2[j] +
	( - x * x + xj * xj + 2 * xj * (x - 2 * xj1) + 2 * xj1 * xj1) * Y2[j+1] ) ) / ( 24.0 * (xj - xj1)) ;
	return value ;
}

template<class T>
void Spline_Interp<T>::copy(const Spline_Interp<T> &spline)
{
	Base_Interp<T>::copy(spline);
	Y2 = spline.Y2 ;
}

template<class T>
Spline_Interp<T>& Spline_Interp<T>::operator= (const Spline_Interp<T> & spline)
{
	// We don't want to waste resources 
	if (this == &spline) {return *this;}
	copy(spline);
	return *this;
}

template<class T>
Spline_Interp<T>::~Spline_Interp() {}

#endif