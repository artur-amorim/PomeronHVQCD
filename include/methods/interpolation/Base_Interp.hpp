#ifndef BASE_INTERP_HPP
#define BASE_INTERP_HPP

#include <vector>
#include <cmath>
#include <stdexcept>

template<class T>
class Base_Interp{
    private:
        int n, m, jsav, cor, dj;
	    std::vector<T> XX, YY;
	protected:
		void copy(const Base_Interp<T> &base);
		int getN() const ;
		int getM() const ;
		int getCor() const;
		T X(const int i) const;
		T Y(const int i) const;
		int locate(const T x);
        int hunt(const T x);
		virtual T rawinterp(int jlo, T x);
        virtual T rawder1(int jl0, T x);
		virtual T rawder2(int jl0, T x);
    public:
        Base_Interp();
        Base_Interp(const std::vector<T> &X, const std::vector<T> &Y, const int mm);
		Base_Interp(const Base_Interp<T> &base);
        T interp(const T x);
        T der1(const T x);
        T der2(const T x);
        Base_Interp<T>& operator= (const Base_Interp<T> &base);
        virtual ~Base_Interp();
};

template<class T>
Base_Interp<T>::Base_Interp():
			n(0), m(0), jsav(0), cor(0), dj(0), XX({}), YY{{}} {}

template<class T>
Base_Interp<T>::Base_Interp(const std::vector<T> &X, const std::vector<T> &Y, const int mm): 
			n(X.size()), m(mm), jsav(0), cor(0), dj(std::min(1,(int)std::pow((T)n,0.25))), XX(X), YY(Y){}


template<class T>
Base_Interp<T>::Base_Interp(const Base_Interp<T> &base):
	n(base.n), m(base.m), jsav(base.jsav), cor(base.cor), dj(base.dj), XX(base.XX), YY(base.YY) {}

template<class T>		
int Base_Interp<T>::locate(const T x)
{
	int ju, jm, jl;
	if ( n < 2 || m < 2 || m > n) throw std::runtime_error("Base_Interp::locate size error");
	bool ascnd = ( XX[n-1] >= XX[0]) ;
	jl = 0;
	ju = n - 1 ;
	while (ju - jl > 1)
	{
		jm = (ju + jl) / 2 ;
		if ( x >= XX[jm] == ascnd) jl = jm ;
		else ju = jm;
	}
	cor = std::abs(jl - jsav) > dj ? 0 : 1 ;
	jsav = jl;
	return std::max(0, std::min(n-m, jl - (m-2)/2)) ;

}

template<class T>
int Base_Interp<T>::hunt(const T x)
{
	int jl = jsav, jm, ju, inc = 1;
	if ( n < 2 || m < 2 || m > n) throw std::runtime_error("Base_Interp::hunt size error");
	bool ascnd = (XX[n-1] >= XX[0]) ;
	// Check if input guess is useful
	if ( jl < 0 || jl > n - 1)
	{
		jl = 0 ;
		ju = n - 1;
	}
	else
	{
		if ( x >= XX[jl] == ascnd )
		{
			for(;;)
			{
				ju = jl + inc ;
				if( ju >= n-1 )
				{
					ju = n - 1 ;
					break ;
				}
				else if ( x < XX[ju] == ascnd ) break ;
				else
				{
					jl = ju ;
					inc += inc ;
				}
			}
		}
		else
		{
			ju = jl;
			for(;;)
			{
				jl = jl - inc ;
				if (jl <= 0)
				{
					jl = 0 ;
					break ;
				}
				else if ( x >= XX[jl] == ascnd) break;
				else
				{
					ju = jl;
					inc += inc;
				}
			}
		}
	}
	while ( ju - jl > 1)
	{
		jm = (ju + jl) / 2 ;
		if( x >= XX[jm] == ascnd) jl = jm ;
		else ju = jm ;
	}
	cor = std::abs( jl -jsav ) > dj ? 0 : 1;
	jsav = jl ;
	return std::max(0, std::min(n-m, jl - (m-2)/2 ));
}

template<class T>
void Base_Interp<T>::copy(const Base_Interp<T> &base)
{
	n    = base.n;
	m    = base.m;
	jsav = base.jsav;
	cor  = base.cor;
	dj   = base.dj;
	XX   = base.XX;
	YY   = base.YY;
}

template<class T>
int Base_Interp<T>::getN() const
{
	return this->n;
}

template<class T>
int Base_Interp<T>::getM() const
{
	return this->m;
}

template<class T>
int Base_Interp<T>::getCor() const
{
	return this->cor;
}

template<class T>
T Base_Interp<T>::interp( const T x )
{
	// Given x, return interpolated value, using data pointed to xx and yy
	int jlo = cor ? hunt(x) : locate(x) ;
	return rawinterp(jlo, x) ;
}

template<class T>
T Base_Interp<T>::der1( const T x )
{
	// Given x, return interpolated value, using data pointed to xx and yy
	int jlo = cor ? hunt(x) : locate(x) ;
	return rawder1(jlo, x) ;
}

template<class T>
T Base_Interp<T>::der2( const T x )
{
	// Given x, return interpolated value, using data pointed to xx and yy
	int jlo = cor ? hunt(x) : locate(x) ;
	return rawder2(jlo, x) ;
}

template<class T>
T Base_Interp<T>::X(const int i) const
{
	return XX[i];
}

template<class T>
T Base_Interp<T>::Y(const int i) const
{
	return YY[i];
}

template<class T>
T Base_Interp<T>::rawinterp(int jlo, T x) {return 0;}

template<class T>
T Base_Interp<T>::rawder1(int jlo, T x) {return 0;} 

template<class T>
T Base_Interp<T>::rawder2(int jlo, T x) {return 0;} 

template<class T>
Base_Interp<T>& Base_Interp<T>::operator= (const Base_Interp<T> &base)
{
	// We don't want to waste resources 
	if (this == &base) {return *this;}
	copy(base);
	return *this;
}

template<class T>
Base_Interp<T>::~Base_Interp() {}

#endif