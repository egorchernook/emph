#ifndef MAIN_CPP_SOLUTION_T_CONCEPT_HPP
#define MAIN_CPP_SOLUTION_T_CONCEPT_HPP

#include <cassert>
#include <cstdlib>
#include <iostream>

namespace Numerical_methods {
    template<typename T>
    concept solution_t_concept = requires(T lhs, T rhs) {
        lhs + rhs;
        lhs - rhs;
        lhs * 2.0;
        2.0 * rhs;
        lhs / 2.0;
        lhs += rhs;
        lhs -= rhs;
        lhs *= 2.0;
        lhs /= 2.0;
        lhs[0];
        lhs.at(0);
        lhs.size();
    };
    template<typename T>
    concept solution_element_t = requires(T lhs, T rhs) {
        lhs + rhs;
        lhs - rhs;
        lhs * 2.0;
        2.0 * rhs;
        lhs / 2.0;
        lhs += rhs;
        lhs -= rhs;
        lhs *= 2.0;
        lhs /= 2.0;
    };

    template<solution_t_concept solution_t>
    struct return_t {
        solution_t solution;
        double error;
    };
}

#endif //MAIN_CPP_SOLUTION_T_CONCEPT_HPP
