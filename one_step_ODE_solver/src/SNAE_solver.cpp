#include "SNAE_solver.h"

SNAE_solver::~SNAE_solver() noexcept
{
    if (repr_ptr != nullptr) {
        delete repr_ptr;
        repr_ptr = nullptr;
    }
}
