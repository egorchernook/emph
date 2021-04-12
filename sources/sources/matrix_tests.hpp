#ifndef MAIN_CPP_MATRIX_TESTS_HPP
#define MAIN_CPP_MATRIX_TESTS_HPP

#include "matrix.hpp"
namespace Numerical_methods {
    class matrix_tests {
    private:
    public:
        matrix_tests();

        static bool creation_tests();

        static bool two_parameters_at_test();

        static bool direct_assignment_test();

        static bool addition_assignment_test();

        static bool subtraction_assignment_test();

        static bool multiplication_with_number_assignment_test();

        static bool division_by_number_assignment_test();

        static bool addition_test();

        static bool subtraction_test();

        static bool multiplication_with_number_test();

        static bool division_by_number_test();

        static bool two_matrix_multiplication_test();

        static bool inverse_Gauss_test();
    };
}

#endif //MAIN_CPP_MATRIX_TESTS_HPP
