#include "SNAE_solver.h"

SNAE_solver::~SNAE_solver() noexcept
{
    delete repr_ptr;
    repr_ptr = nullptr;
}
