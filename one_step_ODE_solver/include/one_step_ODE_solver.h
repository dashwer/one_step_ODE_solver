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

#define ODE_SOLVER Rosenbrock_two_stage_solver
#define NUMBER_TYPE svr::complex_t

namespace svr
{
    const std::size_t start_steps_amount = 100;
    using real_vec_t = std::vector<double>;
    using real_matr_t = std::vector<real_vec_t>;
    using domain_t = std::pair<double, double>;
    using argument_t = real_vec_t;
    using function_t = std::function<double(argument_t)>;
    using system_t = std::vector<function_t>;
    using complex_t = std::complex<double>;
    using complex_vec_t = std::vector<complex_t>;
    using complex_matr_t = std::vector<complex_vec_t>;
}

/// VARIOUS COEFFICIENTS

struct Coefficients_representation
{
    virtual void virtual_method() const = 0;
    virtual ~Coefficients_representation() {};
};

struct Runge_Kutta_method : Coefficients_representation
{
    virtual void virtual_method() const override {};

    svr::real_matr_t A;
    svr::real_vec_t b;
    svr::real_vec_t c;

    virtual ~Runge_Kutta_method() override {};
};

struct Rosenbrock_one_stage_method : Coefficients_representation
{
    virtual void virtual_method() const override {};

    svr::complex_t alpha;

    virtual ~Rosenbrock_one_stage_method() override {};
};

struct Rosenbrock_two_stage_method : Coefficients_representation
{
    virtual void virtual_method() const override {};

    svr::complex_t alpha_1;
    svr::complex_t alpha_2;
    svr::complex_t b1;
    svr::complex_t b2;
    svr::complex_t a_21;
    svr::complex_t c_21;

    virtual ~Rosenbrock_two_stage_method() override {};
};

/// ODE SOLVER

template <typename T>
class One_step_ODE_solver_representation
{
public:
    One_step_ODE_solver_representation(const svr::system_t& F, svr::domain_t domain, const svr::real_vec_t& u0) :
        F_{F}, domain_{domain}, u0_{u0}, SLAE_solver_ptr_{new SLAE_solver<T>},
        SNAE_solver_ptr_{new SNAE_solver}, jacobian_ptr_{new Jacobian}
    {}

    virtual svr::real_matr_t get_solution() const = 0;
    virtual ~One_step_ODE_solver_representation();
protected:
    svr::system_t F_;
    svr::domain_t domain_;
    svr::real_vec_t u0_;

    Coefficients_representation* coeffs_ptr_;
    SLAE_solver<T>* SLAE_solver_ptr_;
    SNAE_solver* SNAE_solver_ptr_;
    Jacobian* jacobian_ptr_;
};


template <typename T>
One_step_ODE_solver_representation<T>::~One_step_ODE_solver_representation()
{
    if (coeffs_ptr_ != nullptr) {
        delete coeffs_ptr_;
        coeffs_ptr_ = nullptr;
    }

    if (SNAE_solver_ptr_ != nullptr) {
        delete SNAE_solver_ptr_;
        SNAE_solver_ptr_ = nullptr;
    }

    if (SLAE_solver_ptr_ != nullptr) {
        delete SLAE_solver_ptr_;
        SLAE_solver_ptr_ = nullptr;
    }

    if (jacobian_ptr_ != nullptr) {
        delete jacobian_ptr_;
        jacobian_ptr_ = nullptr;
    }
}

class Runge_Kutta_solver : public One_step_ODE_solver_representation<double>
{
public:
    Runge_Kutta_solver(const svr::system_t& F, svr::domain_t domain, const svr::real_vec_t& u0);

    virtual svr::real_matr_t get_solution() const override;
    virtual ~Runge_Kutta_solver() override {};
};


class Rosenbrock_one_stage_solver : public One_step_ODE_solver_representation<svr::complex_t>
{
public:
    Rosenbrock_one_stage_solver(const svr::system_t& F, svr::domain_t domain, const svr::real_vec_t& u0);

    virtual svr::real_matr_t get_solution() const override;
    virtual ~Rosenbrock_one_stage_solver() override {};
};


class Rosenbrock_two_stage_solver : public One_step_ODE_solver_representation<svr::complex_t>
{
public:
    Rosenbrock_two_stage_solver(const svr::system_t& F, svr::domain_t domain, const svr::real_vec_t& u0);

    virtual svr::real_matr_t get_solution() const override;
    virtual ~Rosenbrock_two_stage_solver() override {};
};


class One_step_ODE_solver
{
public:
    One_step_ODE_solver(const svr::system_t& F, const svr::domain_t& domain, const svr::real_vec_t& u0) :
        repr_ptr{new ODE_SOLVER{F, domain, u0}}
    {}

    One_step_ODE_solver_representation<NUMBER_TYPE>* repr_ptr;

    ~One_step_ODE_solver();
};

#endif // One_step_ODE_solver_H_INCLUDED
