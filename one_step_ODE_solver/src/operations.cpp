#include "operations.h"

opr::real_vec_t operator+(const opr::real_vec_t& lhs, const opr::real_vec_t& rhs)
{
    opr::real_vec_t result;

    for (std::size_t i = 0; i < rhs.size(); ++i) {
        result.push_back(lhs[i] + rhs[i]);
    }

    return result;
}


opr::real_vec_t operator*(const double lhs, const opr::real_vec_t& rhs)
{
    opr::real_vec_t result;

    for (std::size_t i = 0; i < rhs.size(); ++i) {
        result.push_back(lhs * rhs[i]);
    }

    return result;
}
