#include <iostream>
#include <vector>
#include <utility>

using function_t = int;

int main()
{
    One_step_ODE_solver solver_obj = One_step_ODE_solver::create(F, domain, u0);
    /// F - vector, domain - pair, u0 - function_t

    U = solver_obj.get_solution(); /// U - vector<function_t>
}
