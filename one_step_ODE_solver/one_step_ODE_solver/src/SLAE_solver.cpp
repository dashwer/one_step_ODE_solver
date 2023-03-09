#include "SLAE_solver.h"

SLAE_solver::~SLAE_solver()
{
    delete repr_ptr;
    repr_ptr = nullptr;
}
