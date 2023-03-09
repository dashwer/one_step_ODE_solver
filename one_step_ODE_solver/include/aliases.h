#ifndef ALIASES_H_INCLUDED
#define ALIASES_H_INCLUDED

#include <functional>
#include "functor.h"
#include <vector>
#include <utility>

#define ODE_SOLVER Runge_Kutta_solver

using real_vec_t = std::vector<double>;
using real_matr_t = std::vector<real_vec_t>;
using function_t = std::function<double(real_vec_t)>;
using domain_t = std::pair<double, double>;
using system_t = std::vector<Functor>;

#endif // ALIASES_H_INCLUDED
