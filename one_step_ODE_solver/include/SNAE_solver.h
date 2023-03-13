#ifndef SNAE_SOLVER_H_INCLUDED
#define SNAE_SOLVER_H_INCLUDED

#include <functional>
#include "jacobian.h"
#include <iostream>
#include "operations.h"
#include "SLAE_solver.h"
#include <vector>

#define SNAE_SOLVER Newton_solver

namespace snae
{
    const std::size_t max_iterations_amount = 100;
    const double eps = 1e-6;
    using real_vec_t = std::vector<double>;
    using real_matr_t = std::vector<real_vec_t>;
    using argument_t = real_vec_t;
    using function_t = std::function<double(argument_t)>;
    using system_t = std::vector<function_t>;
}


class SNAE_solver_representation
{
public:
    SNAE_solver_representation(const snae::system_t& f) :
        f_{f}, SLAE_solver_ptr_{new SLAE_solver<double>},
        jacobian_ptr_{new Jacobian}
    {}

    SNAE_solver_representation() :
        SLAE_solver_ptr_{new SLAE_solver<double>},
        jacobian_ptr_{new Jacobian}
    {}

    void initialize(const snae::system_t& f) { f_ = f; };
    virtual snae::real_vec_t get_solution() const = 0;
    virtual ~SNAE_solver_representation() noexcept;

protected:
    snae::system_t f_;
    SLAE_solver<double>* SLAE_solver_ptr_;
    Jacobian* jacobian_ptr_;
};


class FPI_solver : public SNAE_solver_representation
{
public:
    using SNAE_solver_representation::SNAE_solver_representation;

    virtual snae::real_vec_t get_solution() const override;

    virtual ~FPI_solver() noexcept override {};
};


class Newton_solver : public SNAE_solver_representation
{
public:
    using SNAE_solver_representation::SNAE_solver_representation;

    virtual snae::real_vec_t get_solution() const override;

    virtual ~Newton_solver() noexcept override {};
};


class SNAE_solver
{
public:
    SNAE_solver(const snae::system_t& f) :
        repr_ptr{new SNAE_SOLVER{f}}
    {}

    SNAE_solver() :
        repr_ptr{new SNAE_SOLVER}
    {}

    SNAE_solver_representation* repr_ptr;

    ~SNAE_solver() noexcept;
};

#endif // SNAE_SOLVER_H_INCLUDED
