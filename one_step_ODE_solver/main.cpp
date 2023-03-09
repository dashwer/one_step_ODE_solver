#include "one_step_ODE_solver.h"
#include <iostream>
#include <functional>
#include <vector>
#include <set>
#include <utility>

using real_vec_t = std::vector<double>;
using real_matr_t = std::vector<real_vec_t>;
using argument_t = real_vec_t;
using function_t = std::function<double(argument_t)>;
using domain_t = std::pair<double, double>;
using system_t = std::vector<function_t>;

double f1(const real_vec_t& values)
{
    return 0.0;
}

double f2(const real_vec_t& values)
{
    return 0.0;
}

int main()
{
    const system_t F = {f1, f2};
    const domain_t domain{0.0, 1.0};
    const argument_t u0 = {1.0, 2.0};

    const One_step_ODE_solver solver_obj{F, domain, u0};

    const real_matr_t U = solver_obj.repr_ptr->get_solution();
}
