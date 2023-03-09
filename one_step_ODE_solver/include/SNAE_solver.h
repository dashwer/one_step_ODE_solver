#ifndef SNAE_SOLVER_H_INCLUDED
#define SNAE_SOLVER_H_INCLUDED

#include <functional>
#include <iostream>
#include <vector>

#define SNAE_SOLVER Newton_solver

namespace snae
{
    using real_vec_t = std::vector<double>;
    using argument_t = real_vec_t;
    using function_t = std::function<double(argument_t)>;
    using system_t = std::vector<function_t>;
}


class SNAE_solver_representation
{
public:
    SNAE_solver_representation(const snae::system_t& f) :
        f_{f}
    {}

    SNAE_solver_representation() = default;

    void initialize(const snae::system_t& f) { f_ = f; };
    virtual snae::real_vec_t get_solution() const = 0;
    virtual ~SNAE_solver_representation() noexcept {};

protected:
    snae::system_t f_;
};


class FPI_solver : public SNAE_solver_representation
{
public:
    using SNAE_solver_representation::SNAE_solver_representation;

    virtual snae::real_vec_t get_solution() const override { return snae::real_vec_t{0.0}; };

    virtual ~FPI_solver() noexcept override {};
};


class Newton_solver : public SNAE_solver_representation
{
public:
    using SNAE_solver_representation::SNAE_solver_representation;

    virtual snae::real_vec_t get_solution() const override
    {
        snae::real_vec_t vec;
        vec.resize(f_.size(), 0.0);
        return vec;
    };

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

    SNAE_solver_representation* repr_ptr = nullptr;

    ~SNAE_solver() noexcept;
};

#endif // SNAE_SOLVER_H_INCLUDED
