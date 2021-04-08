#ifndef MAIN_CPP_ODE_SOLVER_HPP
#define MAIN_CPP_ODE_SOLVER_HPP

#include <functional>
#include <type_traits>
#include <utility>

#include "solution_t_concept.hpp"
#include "ring_buffer.hpp"

namespace Numerical_methods{

    template<solution_t_concept solution_t>
    class Cauchy_problem_solver{};

    template<solution_t_concept solution_t>
    class ODE_solver{
    private:

    public:
    };
}

#endif //MAIN_CPP_ODE_SOLVER_HPP
