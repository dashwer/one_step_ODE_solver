#ifndef OPERATIONS_H_INCLUDED
#define OPERATIONS_H_INCLUDED

#include <complex>
#include <vector>

template <typename T>
T operator+(const T& lhs, const T& rhs);

template <typename T>
T operator-=(T& lhs, const T& rhs);

template <typename T>
T operator*(const double lhs, const T& rhs);

template <typename T>
std::vector<T> operator*(const double lhs, const std::vector<T>& rhs);

template <typename T>
std::vector<T>& operator-=(std::vector<T>& lhs, const std::vector<T>& rhs);


template <typename T>
T operator+(const T& lhs, const T& rhs)
{
    T result;

    for (std::size_t i = 0; i < rhs.size(); ++i) {
        result.push_back(lhs[i] + rhs[i]);
    }

    return result;
}

template <typename T>
T operator-=(T& lhs, const T& rhs)
{
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        lhs[i] -= rhs[i];
    }

    return lhs;
}

template <typename T>
T operator*(const double lhs, const T& rhs)
{
    T result;

    for (std::size_t i = 0; i < rhs.size(); ++i) {
        result.push_back(lhs * rhs[i]);
    }

    return result;
}

template <typename T>
std::vector<T> operator*(const double lhs, const std::vector<T>& rhs)
{
    std::vector<T> result;

    for (auto& row : rhs) {
        result.push_back(lhs * row);
    }

    return result;
}

template <typename T>
std::vector<T>& operator-=(std::vector<T>& lhs, const std::vector<T>& rhs)
{
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        lhs[i] -= rhs[i];
    }

    return lhs;
}

#endif // OPERATIONS_H_INCLUDED
