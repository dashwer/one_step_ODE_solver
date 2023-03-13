#include "one_step_ODE_solver.h"
#include <iostream>
#include <functional>
#include <vector>
#include <set>
#include <utility>

#include "SLAE_solver.h"

using real_vec_t = std::vector<double>;
using real_matr_t = std::vector<real_vec_t>;
using argument_t = real_vec_t;
using function_t = std::function<double(argument_t)>;
using domain_t = std::pair<double, double>;
using system_t = std::vector<function_t>;
using complex_t = std::complex<double>;
using complex_vec_t = std::vector<complex_t>;
using complex_matr_t = std::vector<complex_vec_t>;

double f1(const real_vec_t& values)
{
    return values[1];
}

double f2(const real_vec_t& values)
{
    return -values[0];
}

int main()
{
    const system_t F = {f1, f2};
    const domain_t domain{0.0, 10.0};
    const argument_t u0 = {1.0, 0.0};

    const One_step_ODE_solver solver_obj{F, domain, u0};

    const real_matr_t U = solver_obj.repr_ptr->get_solution();

    for (std::size_t i = 0; i < 100; ++i) {
        for (std::size_t j = 0; j < F.size(); ++j) {
            std::cout << U[i][j] << " ";
        }
        std::cout << std::endl;
    }
}
