#include "jacobian.h"

Jacobian::~Jacobian()
{
    if (repr_ptr != nullptr) {
        delete repr_ptr;
        repr_ptr = nullptr;
    }
}
