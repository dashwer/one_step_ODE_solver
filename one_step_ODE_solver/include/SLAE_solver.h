#ifndef SLAE_SOLVER_H_INCLUDED
#define SLAE_SOLVER_H_INCLUDED

#include <cmath>
#include <stdexcept>
#include <vector>

#define SLAE_SOLVER Gauss_solver<T>

namespace slae
{
    const double eps = 1e-12;
    using real_vec_t = std::vector<double>;
    using real_matr_t = std::vector<real_vec_t>;
    template <typename T> using vec_t = std::vector<T>;
    template <typename T> using matr_t = std::vector<vec_t<T>>;
}


template <typename T>
class SLAE_solver_representation
{
public:
    SLAE_solver_representation(const slae::matr_t<T>& A, const slae::vec_t<T>& b) :
        A_{A}, b_{b}
    {}

    SLAE_solver_representation() = default;

    void initialize(const slae::matr_t<T>& A, const slae::vec_t<T>& b)
    { A_ = A; b_ = b; };

    virtual slae::vec_t<T> get_solution() const = 0;

    virtual ~SLAE_solver_representation() noexcept {};
protected:
    slae::matr_t<T> A_;
    slae::vec_t<T> b_;
};


template <typename T>
class Gauss_solver : public SLAE_solver_representation<T>
{
    using SLAE_solver_representation<T>::SLAE_solver_representation;

    virtual slae::vec_t<T> get_solution() const override;

    virtual ~Gauss_solver() override {};
};


template <typename T>
class SLAE_solver
{
public:
    SLAE_solver(const slae::matr_t<T>& A, const slae::vec_t<T>& b) :
        repr_ptr{new SLAE_SOLVER{A, b}}
    {}

    SLAE_solver() :
        repr_ptr{new SLAE_SOLVER}
    {}

    SLAE_solver_representation<T>* repr_ptr = nullptr;

    ~SLAE_solver() noexcept;
};


template <typename T>
slae::vec_t<T> Gauss_solver<T>::get_solution() const
{
    int N = this->A_.size();
    slae::vec_t<T> x;
    x.resize(N);

    slae::matr_t<T> augmented_A;
    augmented_A.resize(N);
    for (int i = 0; i < N; ++i) {
        augmented_A[i] = this->A_[i];
        augmented_A[i].push_back(this->b_[i]);
    }

    for (int i = 0; i < N; i++) {

        if (std::abs(augmented_A[i][i]) < slae::eps) {

            bool zero_col = true;

            for (int j = i+1; j < N; j++) {
                if (std::abs(augmented_A[j][i]) >= slae::eps) {
                    zero_col = false;
                    std::swap(augmented_A[j], augmented_A[i]);
                    break;
                }
            }
            if (zero_col) {
                throw std::runtime_error{"Ill-conditioned matrix"};
            };
        }

        for (int j = i+1; j < N; j++) {
            if (std::abs(augmented_A[j][i]) >= slae::eps) {
                T coeff = augmented_A[j][i]/augmented_A[i][i];
                augmented_A[j][i] = 0.0;

                for (int k = i+1; k < N+1; k++) {
                    augmented_A[j][k] -= coeff * augmented_A[i][k];
                }
            }
        }
    }

    for (int i = N-1; i >= 0; --i) {
        T sum = 0.0;
        for (int j = i+1; j < N; ++j) {
            sum += augmented_A[i][j] * x[j];
        }
        x[i] = augmented_A[i][N] - sum;
        x[i] /= augmented_A[i][i];
    }

    return x;
}


template <typename T>
SLAE_solver<T>::~SLAE_solver() noexcept
{
    if (repr_ptr != nullptr) {
        delete repr_ptr;
        repr_ptr = nullptr;
    }
}

#endif // SLAE_SOLVER_H_INCLUDED
