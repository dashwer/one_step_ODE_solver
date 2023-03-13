#ifndef OPERATIONS_H_INCLUDED
#define OPERATIONS_H_INCLUDED

#include <complex>
#include <cmath>
#include <iostream>
#include <vector>

namespace opr
{
    using real_vec_t = std::vector<double>;
    using real_matr_t = std::vector<real_vec_t>;
    using complex_t = std::complex<double>;
    using complex_vec_t = std::vector<complex_t>;
    using complex_matr_t = std::vector<complex_vec_t>;
}

template <typename T>
T operator+(const T& lhs, const T& rhs)
{
    T result;

    for (std::size_t i = 0; i < lhs.size(); ++i) {
        result.push_back(lhs[i] + rhs[i]);
    }

    return result;
}

template <typename T>
T& operator+=(T& lhs, const T& rhs)
{
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        lhs[i] += rhs[i];
    }

    return lhs;
}

template <typename T>
T operator-(const T& lhs, const T& rhs)
{
    T result;

    for (std::size_t i = 0; i < lhs.size(); ++i) {
        result.push_back(lhs[i] - rhs[i]);
    }

    return result;
}

template <typename T>
T& operator-=(T& lhs, const T& rhs)
{
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        lhs[i] -= rhs[i];
    }

    return lhs;
}

template <typename T>
opr::real_vec_t operator*(const double lhs, const std::vector<T>& rhs)
{
    std::vector<double> result;

    for (std::size_t i = 0; i < rhs.size(); ++i) {
        result.push_back(lhs * rhs[i]);
    }

    return result;
}

template <typename T>
opr::complex_vec_t operator*(const opr::complex_t lhs, const std::vector<T>& rhs)
{
    std::vector<opr::complex_t> result;

    for (std::size_t i = 0; i < rhs.size(); ++i) {
        result.push_back(lhs * rhs[i]);
    }

    return result;
}

template <typename T>
std::vector<T> linear_comb(const std::vector<T>& coeffs, const std::vector<T>& vectors)
{
    std::size_t vector_size = vectors.size() / coeffs.size();
    std::vector<T> result;

    for (std::size_t i = 0; i < vector_size; ++i) {
        T elem = 0.0;
        for (std::size_t j = 0; j < coeffs.size(); ++j) {
            elem += coeffs[j] * vectors[j * vector_size + i];
        }
        result.push_back(elem);
    }

    return result;
}

double norm(const opr::real_vec_t& vec);

opr::real_vec_t real(const opr::complex_vec_t& vec);


#endif // OPERATIONS_H_INCLUDED
