#include "jacobian.h"

Jacobian::~Jacobian()
{
    if (repr_ptr != nullptr) {
        delete repr_ptr;
        repr_ptr = nullptr;
    }
}


jcb::real_matr_t Jacobian_simple::get_matrix(const jcb::real_vec_t& values) const
{
    jcb::real_matr_t output(values.size());

    for (std::size_t i = 0; i < output.size(); ++i) {
        output[i].resize(output.size());

        for (std::size_t j = 0; j < output.size(); ++j) {
            jcb::real_vec_t values_up = values;
            values_up[j] += jcb::h;
            jcb::real_vec_t values_down = values;
            values_down[j] -= jcb::h;
            output[i][j] = (F_[i](values_up) - F_[i](values_down)) / (2.0 * jcb::h);
        }
    }

    return output;
}
