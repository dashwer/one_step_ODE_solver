#include "one_step_ODE_solver.h"

Runge_Kutta_solver::Runge_Kutta_solver(const svr::system_t& F, svr::domain_t domain, const svr::real_vec_t& u0) :
    One_step_ODE_solver_representation::One_step_ODE_solver_representation(F, domain, u0)
{
    coeffs_ptr_ = new Runge_Kutta_method;
    Runge_Kutta_method* coeffs_ptr = dynamic_cast<Runge_Kutta_method*>(coeffs_ptr_);

    coeffs_ptr->A = {{5.0/12.0, -1.0/12.0},
                      {3.0/4.0, 1.0/4.0}};

    coeffs_ptr->b = {3.0/4.0, 1.0/4.0};
    coeffs_ptr->c = {1.0/3.0, 1.0};

    coeffs_ptr = nullptr;
}


svr::real_matr_t Runge_Kutta_solver::get_solution() const
{
    Runge_Kutta_method* coeffs_ptr = dynamic_cast<Runge_Kutta_method*>(coeffs_ptr_);

    const double tau = (domain_.second - domain_.first) / svr::start_steps_amount;

    svr::real_matr_t U;
    U.resize(svr::start_steps_amount);
    U[0] = u0_;

    const std::size_t s = coeffs_ptr->b.size();
    const std::size_t N = F_.size();

    for (std::size_t n = 0; n < svr::start_steps_amount; ++n) {
        const double t_n = domain_.first + tau * n;
        svr::system_t f;

        for (std::size_t i = 0; i < s; ++i) {
            for (std::size_t j = 0; j < N; ++j) {

                auto f_elem = [=, &coeffs_ptr](const svr::real_vec_t& k)
                {
                    svr::real_vec_t args = {t_n + coeffs_ptr->c[i] * tau};
                    svr::real_vec_t a_k = linear_comb(coeffs_ptr->A[i], k);

                    for (std::size_t p = 0; p < N; ++p) {
                        const double value = U[n][p] + tau * a_k[p];
                        args.push_back(value);
                    }

                    return F_[j](args) - k[i*s + j];
                };

                f.emplace_back(f_elem);
            }
        }

        SNAE_solver_ptr_->repr_ptr->initialize(f);
        svr::real_vec_t k = SNAE_solver_ptr_->repr_ptr->get_solution();

        svr::real_vec_t b_k = linear_comb(coeffs_ptr->b, k);
        U[n+1] = U[n] + tau * b_k;
    }
    coeffs_ptr = nullptr;

    return U;
}


Rosenbrock_one_stage_solver::Rosenbrock_one_stage_solver(const svr::system_t& F, svr::domain_t domain, const svr::real_vec_t& u0) :
    One_step_ODE_solver_representation::One_step_ODE_solver_representation(F, domain, u0)
{
    coeffs_ptr_ = new Rosenbrock_one_stage_method;
    Rosenbrock_one_stage_method* coeffs_ptr = dynamic_cast<Rosenbrock_one_stage_method*>(coeffs_ptr_);

    coeffs_ptr->alpha = svr::complex_t(0.5, 0.5);
}

svr::real_matr_t Rosenbrock_one_stage_solver::get_solution() const
{
    Rosenbrock_one_stage_method* coeffs_ptr = dynamic_cast<Rosenbrock_one_stage_method*>(coeffs_ptr_);

    const double tau = (domain_.second - domain_.first) / svr::start_steps_amount;

    svr::real_matr_t U;
    U.resize(svr::start_steps_amount);
    U[0] = u0_;
    const std::size_t N = u0_.size();
    jacobian_ptr_->repr_ptr->initialize(F_);

    for (std::size_t n = 0; n < svr::start_steps_amount; ++n) {
        svr::real_matr_t J = jacobian_ptr_->repr_ptr->get_matrix(U[n]);
        svr::complex_matr_t A;
        A.resize(N);
        svr::complex_vec_t b;

        for (std::size_t i = 0; i < N; ++i) {
            A[i].resize(N, svr::complex_t(0.0, 0.0));
            A[i][i] = 1.0;
            A[i] -= coeffs_ptr->alpha * tau * J[i];
            b.push_back(F_[i](U[n]));
        }

        SLAE_solver_ptr_->repr_ptr->initialize(A, b);
        svr::complex_vec_t k = SLAE_solver_ptr_->repr_ptr->get_solution();

        U[n+1] = U[n] + tau * real(k);
    }

    return U;
}


Rosenbrock_two_stage_solver::Rosenbrock_two_stage_solver(const svr::system_t& F, svr::domain_t domain, const svr::real_vec_t& u0) :
    One_step_ODE_solver_representation::One_step_ODE_solver_representation(F, domain, u0)
{
    coeffs_ptr_ = new Rosenbrock_two_stage_method;
    Rosenbrock_two_stage_method* coeffs_ptr = dynamic_cast<Rosenbrock_two_stage_method*>(coeffs_ptr_);

    coeffs_ptr->alpha_1 = svr::complex_t(1.0/6.0, 1.0/6.0);
    coeffs_ptr->alpha_2 = 0.25;
    coeffs_ptr->b1 = svr::complex_t(23.0/39.0, -2.0/39.0);
    coeffs_ptr->b2 = 16.0/39.0;
    coeffs_ptr->a_21 = svr::complex_t(0.5, -0.875);
    coeffs_ptr->c_21 = svr::complex_t(0.75, -11.0/48.0);
}

svr::real_matr_t Rosenbrock_two_stage_solver::get_solution() const
{
    Rosenbrock_two_stage_method* coeffs_ptr = dynamic_cast<Rosenbrock_two_stage_method*>(coeffs_ptr_);

    const double tau = (domain_.second - domain_.first) / svr::start_steps_amount;

    svr::real_matr_t U;
    U.resize(svr::start_steps_amount);
    U[0] = u0_;
    const std::size_t N = u0_.size();
    jacobian_ptr_->repr_ptr->initialize(F_);

    for (std::size_t n = 0; n < svr::start_steps_amount; ++n) {
        svr::real_matr_t J = jacobian_ptr_->repr_ptr->get_matrix(U[n]);
        svr::complex_matr_t A;
        A.resize(N);
        svr::complex_vec_t b;

        for (std::size_t i = 0; i < N; ++i) {
            A[i].resize(N, 0.0);
            A[i][i] = 1.0;
            A[i] -= coeffs_ptr->alpha_1 * tau * J[i];
            b.push_back(F_[i](U[n]));
        }

        SLAE_solver_ptr_->repr_ptr->initialize(A, b);
        svr::complex_vec_t k1 = SLAE_solver_ptr_->repr_ptr->get_solution();

        J = jacobian_ptr_->repr_ptr->get_matrix(U[n] + tau * real(coeffs_ptr->a_21 * k1));
        A.resize(N);
        b.clear();

        for (std::size_t i = 0; i < N; ++i) {
            A[i].resize(N, 0.0);
            A[i][i] = 1.0;
            A[i] -= coeffs_ptr->alpha_2 * tau * J[i];
            b.push_back(F_[i](U[n] + tau * real(coeffs_ptr->c_21 * k1)));
        }

        SLAE_solver_ptr_->repr_ptr->initialize(A, b);
        svr::complex_vec_t k2 = SLAE_solver_ptr_->repr_ptr->get_solution();

        U[n+1] = U[n] + tau * real(coeffs_ptr->b1 * k1 + coeffs_ptr->b2 * k2);
    }

    return U;
}


One_step_ODE_solver::~One_step_ODE_solver()
{
    if (repr_ptr != nullptr) {
        delete repr_ptr;
        repr_ptr = nullptr;
    }
}
