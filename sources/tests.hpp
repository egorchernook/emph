#ifndef MAIN_CPP_TESTS_HPP
#define MAIN_CPP_TESTS_HPP

#include "SNAE_solver_tests.hpp"
#include "matrix_tests.hpp"
#include "vector_tests.hpp"
#include "ring_buffer_tests.hpp"
#include "one_step_methods_tests.hpp"
#include "multisteps_methods_tests.hpp"

namespace Numerical_methods {
    inline void test() {
        Numerical_methods::vector_tests();
        Numerical_methods::matrix_tests();
        Numerical_methods::SNAE_solver_tests();
        Numerical_methods::ring_buffer_tests();
        Numerical_methods::one_step_methods_tests();
        Numerical_methods::multisteps_methods_tests();
    }
}
#endif //MAIN_CPP_TESTS_HPP
