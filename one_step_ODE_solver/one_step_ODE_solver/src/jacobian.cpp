#include "jacobian.h"

Jacobian::~Jacobian()
{
    delete repr_ptr;
    repr_ptr = nullptr;
}
