#ifndef ONE_STEP_ODE_SOLVER_H_INCLUDED
#define ONE_STEP_ODE_SOLVER_H_INCLUDED

#include <vector>
#include <utility>

using domain_t = std::pair<double, double>;
using system_t = std::vector<functor_t>;
using real_vec_t = std::vector<double>;
using real_matr_t = std::vector<real_vec_t>;

template <typename T>
class one_step_ODE_solver
{
public:
    using output_t = std::vector<T>;

    one_step_ODE_solver() = default;
    static one_step_ODE_solver<T> create(const system_t& F, const domain_t& domain, const T& u0);

    output_t get_solution();

    ~one_step_ODE_solver() noexcept;

private:
    one_step_ODE_solver_repr* repr_ptr;
};


class one_step_ODE_solver_repr
{
public:
    one_step_ODE_solver_repr(const system_t& F, domain_t domain, const real_vec_t& u0);

    virtual real_matr_t compute_solution() = 0;
protected:
    system_t F_;
    domain_t domain_;
    real_vec_t u0_;
};

class Runge_Kutta_solver : public one_step_ODE_solver_repr
{
public:
    using one_step_ODE_solver_repr::one_step_ODE_solver_repr;

    virtual real_matr_t compute_solution() override;
};


class Rosenbrock_solver : public one_step_ODE_solver_repr
{
public:
    using one_step_ODE_solver_repr::one_step_ODE_solver_repr;

    virtual real_matr_t compute_solution() override;
};

#endif // ONE_STEP_ODE_SOLVER_H_INCLUDED
