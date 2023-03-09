#ifndef One_step_ODE_solver_H_INCLUDED
#define One_step_ODE_solver_H_INCLUDED

#include <complex>
#include <functional>
#include <iostream>
#include "jacobian.h"
#include "operations.h"
#include <vector>
#include "SNAE_solver.h"
#include "SLAE_solver.h"
#include <utility>

#define ODE_SOLVER Runge_Kutta_solver

namespace svr
{
    const std::size_t start_steps_amount = 100;
    using real_vec_t = std::vector<double>;
    using real_matr_t = std::vector<real_vec_t>;
    using domain_t = std::pair<double, double>;
    using argument_t = real_vec_t;
    using function_t = std::function<double(argument_t)>;
    using system_t = std::vector<function_t>;
    using complex_vec_t = std::vector<std::complex<double>>;
}

struct Coefficients_representation
{
    virtual void virtual_method() const = 0;
    virtual ~Coefficients_representation() noexcept {};
};

struct Runge_Kutta_method : Coefficients_representation
{
    virtual void virtual_method() const override {};

    svr::real_matr_t A;
    svr::real_vec_t b;
    svr::real_vec_t c;

    virtual ~Runge_Kutta_method() noexcept override {};
};

struct Rosenbrock_one_stage_method : Coefficients_representation
{
    virtual void virtual_method() const override {};

    double alpha;

    virtual ~Rosenbrock_one_stage_method() noexcept override {};
};

struct Rosenbrock_two_stage_method : Coefficients_representation
{
    virtual void virtual_method() const override {};

    double alpha_1;
    double alpha_2;
    double b1;
    double b2;
    double a_21;
    double c_21;

    virtual ~Rosenbrock_two_stage_method() noexcept override {};
};


class One_step_ODE_solver_repr
{
public:
    One_step_ODE_solver_repr(const svr::system_t& F, svr::domain_t domain, const svr::real_vec_t& u0) :
        F_{F}, domain_{domain}, u0_{u0}, jacobian_ptr_{new Jacobian{F}}
    {}

    virtual svr::real_matr_t get_solution() = 0;
    virtual ~One_step_ODE_solver_repr() noexcept;
protected:
    svr::system_t F_;
    svr::domain_t domain_;
    svr::real_vec_t u0_;

    Coefficients_representation* coeffs_ptr_ = nullptr;
    SLAE_solver* SLAE_solver_ptr_ = nullptr;
    SNAE_solver* SNAE_solver_ptr_ = nullptr;
    Jacobian* jacobian_ptr_ = nullptr;
};

class Runge_Kutta_solver : public One_step_ODE_solver_repr
{
public:
    Runge_Kutta_solver(const svr::system_t& F, svr::domain_t domain, const svr::real_vec_t& u0);

    virtual svr::real_matr_t get_solution() override;
    virtual ~Runge_Kutta_solver() noexcept override {};
};


class Rosenbrock_one_stage_solver : public One_step_ODE_solver_repr
{
public:
    Rosenbrock_one_stage_solver(const svr::system_t& F, svr::domain_t domain, const svr::real_vec_t& u0);

    virtual svr::real_matr_t get_solution() override;
    virtual ~Rosenbrock_one_stage_solver() noexcept override {};
};


class Rosenbrock_two_stage_solver : public One_step_ODE_solver_repr
{
public:
    Rosenbrock_two_stage_solver(const svr::system_t& F, svr::domain_t domain, const svr::real_vec_t& u0);

    virtual svr::real_matr_t get_solution() override;
    virtual ~Rosenbrock_two_stage_solver() noexcept override {};
};


class One_step_ODE_solver
{
public:
    One_step_ODE_solver(const svr::system_t& F, const svr::domain_t& domain, const svr::real_vec_t& u0) :
        repr_ptr{new ODE_SOLVER{F, domain, u0}}
    {}

    One_step_ODE_solver_repr* repr_ptr;

    ~One_step_ODE_solver() noexcept;
};

#endif // One_step_ODE_solver_H_INCLUDED
