#include "SLAE_solver.h"

SLAE_solver::~SLAE_solver()
{
    if (repr_ptr != nullptr) {
        delete repr_ptr;
        repr_ptr = nullptr;
    }
}
