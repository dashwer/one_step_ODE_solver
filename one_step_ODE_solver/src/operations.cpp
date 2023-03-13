#include "operations.h"

double norm(const opr::real_vec_t& vec)
{
    double result = 0.0;
    for (const auto& elem : vec) {
        if (result < fabs(elem)) {
            result = fabs(elem);
        }
    }

    return result;
}

opr::real_vec_t real(const opr::complex_vec_t& vec)
{
    opr::real_vec_t result;

    for (const auto& elem : vec) {
        result.push_back(elem.real());
    }

    return result;
}
