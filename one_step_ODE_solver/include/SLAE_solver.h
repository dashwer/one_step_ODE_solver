#ifndef SLAE_SOLVER_H_INCLUDED
#define SLAE_SOLVER_H_INCLUDED

#include <vector>

#define SLAE_SOLVER Gauss_solver

namespace slae
{
    using real_vec_t = std::vector<double>;
    using real_matr_t = std::vector<real_vec_t>;
}


class SLAE_solver_representation
{
public:
    SLAE_solver_representation(const slae::real_matr_t& A, const slae::real_vec_t& b) :
        A_{A}, b_{b}
    {}

    SLAE_solver_representation() = default;

    void initialize(const slae::real_matr_t& A, const slae::real_vec_t& b)
    { A_ = A; b_ = b; };

    virtual slae::real_vec_t get_solution() const = 0;

    virtual ~SLAE_solver_representation() noexcept {};
protected:
    slae::real_matr_t A_;
    slae::real_vec_t b_;
};


class Gauss_solver : public SLAE_solver_representation
{
    using SLAE_solver_representation::SLAE_solver_representation;

    virtual slae::real_vec_t get_solution() const override {return b_;};

    virtual ~Gauss_solver() override {};
};


class SLAE_solver
{
public:
    SLAE_solver(const slae::real_matr_t& A, const slae::real_vec_t& b) :
        repr_ptr{new SLAE_SOLVER{A, b}}
    {}

    SLAE_solver() :
        repr_ptr{new SLAE_SOLVER}
    {}

    SLAE_solver_representation* repr_ptr = nullptr;

    ~SLAE_solver() noexcept;
};

#endif // SLAE_SOLVER_H_INCLUDED
