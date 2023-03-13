#include "SNAE_solver.h"


SNAE_solver_representation::~SNAE_solver_representation() noexcept
{
    if (SLAE_solver_ptr_ != nullptr) {
        delete SLAE_solver_ptr_;
        SLAE_solver_ptr_ = nullptr;
    }

    if (jacobian_ptr_ != nullptr) {
        delete jacobian_ptr_;
        jacobian_ptr_ = nullptr;
    }
}


snae::real_vec_t FPI_solver::get_solution() const
{
    snae::real_vec_t x, prev_x, f;
    x.resize(f_.size(), 0.0);
    f = x;

    do {
        prev_x = x;

        for (std::size_t i = 0; i < f_.size(); ++i) {
            f[i] = f_[i](prev_x);
        }

        x = prev_x + f;
    } while (norm(x - prev_x) >= snae::eps);

    return x;
}


snae::real_vec_t Newton_solver::get_solution() const
{
    snae::real_vec_t x, prev_x;
    x.resize(f_.size(), 0.0);
    jacobian_ptr_->repr_ptr->initialize(f_);
    snae::real_matr_t A;
    snae::real_vec_t b;
    b.resize(f_.size());

    do {
        prev_x = x;
        A = jacobian_ptr_->repr_ptr->get_matrix(prev_x);

        for (std::size_t i = 0; i < f_.size(); ++i) {
            b[i] = -f_[i](prev_x);
        }

        SLAE_solver_ptr_->repr_ptr->initialize(A, b);
        x = SLAE_solver_ptr_->repr_ptr->get_solution();
        x += prev_x;
    } while (norm(x - prev_x) >= snae::eps);

    return x;
}


SNAE_solver::~SNAE_solver() noexcept
{
    if (repr_ptr != nullptr) {
        delete repr_ptr;
        repr_ptr = nullptr;
    }
}
