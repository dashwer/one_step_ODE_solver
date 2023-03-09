#ifndef JACOBIAN_H_INCLUDED
#define JACOBIAN_H_INCLUDED

#include <functional>
#include <vector>

#define JACOBIAN Jacobian_simple

namespace jcb
{
    using argument_t = std::vector<double>;
    using function_t = std::function<double(argument_t)>;
    using system_t = std::vector<function_t>;
    using real_vec_t = std::vector<double>;
    using real_matr_t = std::vector<real_vec_t>;
}


class Jacobian_representation
{
public:
    Jacobian_representation(const jcb::system_t& F) :
        F_{F}
    {}

    virtual jcb::real_matr_t get_matrix(const jcb::real_vec_t& values) const = 0;
    virtual ~Jacobian_representation() noexcept {};
protected:
    jcb::system_t F_;
};


class Jacobian_simple : public Jacobian_representation
{
    using Jacobian_representation::Jacobian_representation;

    virtual jcb::real_matr_t get_matrix(const jcb::real_vec_t& values) const override
    {
        std::size_t N = F_.size();
        jcb::real_matr_t output;
        jcb::real_vec_t vec;
        vec.resize(N, 0.0);
        output.resize(N, vec);
        return output;
    };

    virtual ~Jacobian_simple() noexcept override {};
};


class Jacobian
{
public:
    Jacobian(const jcb::system_t& F) :
        repr_ptr{new JACOBIAN{F}}
    {}

    Jacobian_representation* repr_ptr = nullptr;
    ~Jacobian() noexcept;
};

#endif // JACOBIAN_H_INCLUDED
