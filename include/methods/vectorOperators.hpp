#ifndef VECTOROPERATORS_HPP
#define VECTOROPERATORS_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

template<class T>
T sum(const std::vector<T> &v)
{
    // Sum all the elements of the vector v
    T s = 0;
    for(int i = 0; i < v.size(); i++) s += v[i];
    return s;
}

template<class T>
std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b)
{
    // Perform the sum of two vectors
    // First check they have the same size
    if( a.size() != b.size() ) 
    {
        std::cout << "Different size vectors!" << std::endl;
        throw "error" ;
    }
    const int l = a.size() ;
    std::vector<T> sol ;
    for(int i = 0; i < l; i++)
    {
        sol.push_back(a[i] + b[i]) ;
    }
    return sol ;
}

template<class T>
std::vector<T> operator+(const T a, const std::vector<T> &x)
{
    // Perform the sum of a number with a vector a + x
    // First check they have the same size
    if( x.size() == 0 )
    {
        std::cout << "Vector with zero size." << std::endl;
        throw "error" ;
    }
    const int n = x.size() ;
    std::vector<T> ans ;
    for(int i = 0; i < n; i++)
    {
        ans.push_back(a + x[i]) ;
    }
    return ans ;
}

template<class T>
std::vector<T> operator-(const std::vector<T> &a, const std::vector<T> &b)
{
    // Perform the difference of two vectors a - b
    // First check they have the same size
    if( a.size() != b.size() )
    {
        std::cout << "Different size vectors!" << std::endl;
        throw "error" ;
    }
    const int l = a.size() ;
    std::vector<T> sol ;
    for(int i = 0; i < l; i++)
    {
        sol.push_back(a[i] - b[i]) ;
    }
    return sol ;
}

template<class T>
std::vector<T> operator-(const T a, const std::vector<T> &x)
{
    // Performs a - x where a is a number and x a vector
    // First check they have the same size
    if( x.size() == 0 )
    {
        std::cout << "Vector with zero size." << std::endl;
        throw "error" ;
    }
    const int n = x.size() ;
    std::vector<T> ans ;
    for(int i = 0; i < n; i++)
    {
        ans.push_back(a - x[i]) ;
    }
    return ans ;
}

template<class T>
std::vector<T> operator*(const std::vector<T> &a, const std::vector<T> &b)
{
    // Perform the product of the entries of two vectors
    // First we check that they have the same size
    if (a.size() != b.size())
    {
        std::cout << "Factors must have the same size" << std::endl ;
        throw "error" ;
    }
    const int n = a.size() ;
    std::vector<T> ans(n) ;
    for (int i = 0; i < n; i++)
    {
        ans[i] = a[i] * b[i] ;
    }
    return ans;
}

template<class T>
std::vector<T> operator*( const T a, const std::vector<T> &x)
{
    // Perform the product of the entries of two vectors
    // First we check that they have the same size
    if ( x.size() == 0) 
    {
        std::cout << "Null vector" << std::endl;
        throw "error" ;
    }
    const int n = x.size() ;
    std::vector<T> ans(n) ;
    for (int i = 0; i < n; i++)
    {
        ans[i] = a * x[i] ;
    }
    return ans;
}

template<class T>
std::vector<T> operator/(const std::vector<T> &num, const std::vector<T> &denom)
{
    // Perform the division of two vectors
    // First we check that they have the same size
    if (num.size() != denom.size())
    {
        std::cout << "Numerator and denominator must have the same size" << std::endl ;
        throw "error" ;
    }
    const int n = num.size() ;
    std::vector<T> ans(n) ;
    for (int i = 0; i < n; i++)
    {
        if (denom[i] == 0)
        {
            std::cout << "Division by 0 found!" << std::endl;
            throw "error";
        }
        else ans[i] = num[i] / denom[i] ;
    }
    return ans;
}

template<class T>
std::vector<T> operator/(const T a , const std::vector<T> &denom)
{
    // Perform the division of a number by a vector
    // First we check that they have the same size
    if (denom.size() == 0) 
    {
        std::cout << "Denominator has 0 size" << std::endl;
        throw "error";
    }
    const int n = denom.size() ;
    std::vector<T> ans(n) ;
    for (int i = 0; i < n; i++)
    {
        if (denom[i] == 0) throw("Division by 0 found!") ;
        else ans[i] = a / denom[i] ;
    }
    return ans;
}

template<class T>
std::vector<T> operator/(const std::vector<T> &num, const T denom)
{
    // Perform the division of a number by a vector
    // First we check that they have the same size
    if (num.size() == 0) throw std::runtime_error("0 size std::vector found in operator/.");
    if (denom == 0) throw std::runtime_error("Dividing std::vector by zero");
    const int n = num.size() ;
    std::vector<T> ans(n) ;
    for (int i = 0; i < n; i++) ans[i] = num[i] / denom ;
    return ans;
}

template<class T>
T norm(const std::vector<T> &x)
{
    // Computes the norm of the vector
    if (x.size() == 0)
    {
        std::cout << "Null vector." << std::endl;
        throw "error";
    }
    T sol = 0;
    const int n = x.size() ;
    for(int i = 0; i < n; i++)
    {
        sol += x[i] * x[i] ;
    }
    return std::sqrt(sol) ;
}

template<class T>
std::vector<T> abs(const std::vector<T> &x)
{
    const int n = x.size();
    std::vector<T> ans(n);
    for (int i = 0; i < n; i++) ans[i] = std::fabs(x[i]);
    return ans;
}

template<class T>
std::vector<T> maximum(const std::vector<T> &x, const std::vector<T> &y)
{
    // Check the size of x is equal to the size of y
    if (x.size() != y.size()) throw std::runtime_error("x and y arguments in maximum don't have the same size.");
    const int n = x.size();
    std::vector<T> ans(n);
    for (int i = 0; i < n; i++) ans[i] = std::max(x[i], y[i]);
    return ans;
}

template<class T>
std::vector<T> sqrt(const std::vector<T> &x)
{
    // Computes the log of a vector
    if (x.size() == 0)
    {
        std::cout << "No element to compute sqrt." << std::endl;
        throw "error";
    }
    const int n = x.size() ;
    std::vector<T> ans(n) ;
    for(int i = 0; i < n; i++)
    {
        if( x[i] < 0)
        {
            std::cout << "Square root of negative number found." << std::endl;
            throw "error";
        }
        else ans[i] = std::sqrt(x[i]) ;
    }
    return ans ;
}

template<class T>
std::vector<T> log(const std::vector<T> &x)
{
    // Computes the log of a vector
    if (x.size() == 0) 
    {
        std::cout << "No element to compute log." << std::endl;
        throw "error";
    }
    const int n = x.size() ;
    std::vector<T> ans(n) ;
    for(int i = 0; i < n; i++)
    {
        if( x[i] <= 0)
        {
            std::cout << "log of negative number or zero found." << std::endl;
            throw "error";
        }
        else ans[i] = std::log(x[i]) ;
    }
    return ans ;
}

template<class T>
std::vector<T> exp(const std::vector<T> &x)
{
    // Computes the exponential of a vector
    if (x.size() == 0) 
    {
        std::cout << "No element to compute exp." << std::endl ;
        throw "error";
    }
    const int n = x.size() ;
    std::vector<T> ans(n) ;
    for(int i = 0; i < n; i++)
    {
        ans[i] = std::exp(x[i]) ;
    }
    return ans ;
}
;
#endif