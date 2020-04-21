#ifndef SEARCH_HPP
#define SEARCH_HPP

#include <iostream>
#include <cmath>
#include <vector>

template<class T>
T binary_search(const std::vector<T> &A, const T a)
{
    int L = 0;
    int R = A.size() - 1;
    while(L <= R)
    {
        int M = std::floor((L + R) / 2.0);
        if (A[M] < a) L = M + 1 ;
        else if (A[M] > a) R = M - 1;
        else return A[M];
    }
    std::cout << "Binary search algorithm failed to find object." << std::endl;
    throw "error";
}

;
#endif