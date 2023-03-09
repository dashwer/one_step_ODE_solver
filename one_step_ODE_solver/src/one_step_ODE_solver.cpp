#include "one_step_ODE_solver.h"

One_step_ODE_solver::~One_step_ODE_solver() noexcept
{
    delete repr_ptr;
    repr_ptr = nullptr;
}


One_step_ODE_solver_repr::~One_step_ODE_solver_repr()
{
    delete coeffs_ptr_;
    coeffs_ptr_ = nullptr;

    delete jacobian_ptr_;
    jacobian_ptr_ = nullptr;
}


Runge_Kutta_solver::Runge_Kutta_solver(const svr::system_t& F, svr::domain_t domain, const svr::real_vec_t& u0) :
    One_step_ODE_solver_repr::One_step_ODE_solver_repr(F, domain, u0)
{
    coeffs_ptr_ = new Runge_Kutta_method;
    Runge_Kutta_method* method_ptr = dynamic_cast<Runge_Kutta_method*>(coeffs_ptr_);

    method_ptr->A = {{5.0/12.0, -1.0/12.0},
                      {3.0/4.0, 1.0/4.0}};

    method_ptr->b = {3.0/4.0, 1.0/4.0};
    method_ptr->c = {1.0/3.0, 1.0};
}


svr::real_matr_t Runge_Kutta_solver::get_solution()
{
    Runge_Kutta_method* method_ptr = dynamic_cast<Runge_Kutta_method*>(coeffs_ptr_);

    const double tau = (domain_.second - domain_.first) / svr::start_steps_amount;

    svr::real_matr_t U;
    U.resize(svr::start_steps_amount);
    U[0] = u0_;

    const std::size_t s = method_ptr->b.size();
    const std::size_t N = F_.size();

    auto linear_comb = [](const svr::real_vec_t& coeffs, const svr::real_vec_t& vectors)
    {
        std::size_t vector_size = vectors.size() / coeffs.size();
        svr::real_vec_t result;

        for (std::size_t i = 0; i < vector_size; ++i) {
            double elem = 0.0;
            for (std::size_t j = 0; j < coeffs.size(); ++j) {
                elem += coeffs[j] * vectors[j * vector_size + i];
            }
            result.push_back(elem);
        }

        return result;
    };

    for (std::size_t n = 0; n < svr::start_steps_amount; ++n) {
        const double t_n = domain_.first + tau * n;
        svr::system_t f;

        for (std::size_t i = 0; i < s; ++i) {
            for (std::size_t j = 0; j < N; ++j) {

                auto f_elem = [&](const svr::real_vec_t& k)
                {
                    svr::real_vec_t args = {t_n + method_ptr->c[i] * tau};
                    svr::real_vec_t a_k = linear_comb(method_ptr->A[i], k);

                    for (std::size_t p = 0; p < N; ++p) {
                        const double value = U[n][p] + tau * a_k[p];
                        args.push_back(value);
                    }

                    double result = F_[j](args) - k[i*s + j];
                    return result;
                };

                f.push_back(f_elem);
            }
        }

        SNAE_solver_ptr_ = new SNAE_solver{f};
        svr::real_vec_t k = SNAE_solver_ptr_->repr_ptr->get_solution();

        if (n+1 != svr::start_steps_amount) {
            delete SNAE_solver_ptr_;
        }

        svr::real_vec_t b_k = linear_comb(method_ptr->b, k);
        U[n+1] = U[n] + tau * b_k;
    }

    method_ptr = nullptr;

    return U;
}


Rosenbrock_one_stage_solver::Rosenbrock_one_stage_solver(const svr::system_t& F, svr::domain_t domain, const svr::real_vec_t& u0) :
    One_step_ODE_solver_repr::One_step_ODE_solver_repr(F, domain, u0)
{
    coeffs_ptr_ = new Rosenbrock_one_stage_method;
    Rosenbrock_one_stage_method* method_ptr = dynamic_cast<Rosenbrock_one_stage_method*>(coeffs_ptr_);

    method_ptr->alpha = 0.5;
}

svr::real_matr_t Rosenbrock_one_stage_solver::get_solution()
{
    Rosenbrock_one_stage_method* method_ptr = dynamic_cast<Rosenbrock_one_stage_method*>(coeffs_ptr_);

    const double tau = (domain_.second - domain_.first) / svr::start_steps_amount;

    svr::real_matr_t U;
    U.resize(svr::start_steps_amount);
    U[0] = u0_;
    const std::size_t N = u0_.size();

    for (std::size_t n = 0; n < svr::start_steps_amount; ++n) {
        svr::real_matr_t J = jacobian_ptr_->repr_ptr->get_matrix(U[n]);
        svr::real_matr_t A;
        A.resize(N);
        svr::real_vec_t b;

        for (std::size_t i = 0; i < N; ++i) {
            A[i].resize(N, 0.0);
            A[i][i] = 1.0;
            A[i] -= method_ptr->alpha * tau * J[i];
            b.push_back(F_[i](U[n]));
        }

        SLAE_solver_ptr_ = new SLAE_solver{A, b};
        svr::real_vec_t k = SLAE_solver_ptr_->repr_ptr->get_solution();

        U[n+1] = U[n] + tau * k;
    }

    return U;
}


Rosenbrock_two_stage_solver::Rosenbrock_two_stage_solver(const svr::system_t& F, svr::domain_t domain, const svr::real_vec_t& u0) :
    One_step_ODE_solver_repr::One_step_ODE_solver_repr(F, domain, u0)
{
    coeffs_ptr_ = new Rosenbrock_two_stage_method;
    Rosenbrock_two_stage_method* method_ptr = dynamic_cast<Rosenbrock_two_stage_method*>(coeffs_ptr_);

    method_ptr->alpha_1 = 0.5;
    method_ptr->alpha_2 = 0.5;
    method_ptr->b1 = 0.5;
    method_ptr->b2 = 0.5;
    method_ptr->a_21 = 0.5;
    method_ptr->c_21 = 0.5;
}

svr::real_matr_t Rosenbrock_two_stage_solver::get_solution()
{
    Rosenbrock_two_stage_method* method_ptr = dynamic_cast<Rosenbrock_two_stage_method*>(coeffs_ptr_);

    const double tau = (domain_.second - domain_.first) / svr::start_steps_amount;

    svr::real_matr_t U;
    U.resize(svr::start_steps_amount);
    U[0] = u0_;
    const std::size_t N = u0_.size();

    for (std::size_t n = 0; n < svr::start_steps_amount; ++n) {
        svr::real_matr_t J = jacobian_ptr_->repr_ptr->get_matrix(U[n]);
        svr::real_matr_t A;
        A.resize(N);
        svr::real_vec_t b;

        for (std::size_t i = 0; i < N; ++i) {
            A[i].resize(N, 0.0);
            A[i][i] = 1.0;
            A[i] -= method_ptr->alpha_1 * tau * J[i];
            b.push_back(F_[i](U[n]));
        }

        SLAE_solver_ptr_ = new SLAE_solver{A, b};
        svr::real_vec_t k1 = SLAE_solver_ptr_->repr_ptr->get_solution();
        delete SLAE_solver_ptr_;

        J = jacobian_ptr_->repr_ptr->get_matrix(U[n] + tau * method_ptr->a_21 * k1);
        A.resize(N);
        b.clear();

        for (std::size_t i = 0; i < N; ++i) {
            A[i].resize(N, 0.0);
            A[i][i] = 1.0;
            A[i] -= method_ptr->alpha_2 * tau * J[i];
            b.push_back(F_[i](U[n] + tau * method_ptr->c_21 * k1));
        }

        SLAE_solver_ptr_ = new SLAE_solver{A, b};
        svr::real_vec_t k2 = SLAE_solver_ptr_->repr_ptr->get_solution();
        delete SLAE_solver_ptr_;

        U[n+1] = U[n] + tau * (method_ptr->b1 * k1 + method_ptr->b2 * k2);
    }

    return U;
}
