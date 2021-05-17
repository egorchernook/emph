#ifndef MAIN_CPP_VECTOR_TESTS_HPP
#define MAIN_CPP_VECTOR_TESTS_HPP

#include "vector.hpp"
namespace Numerical_methods {
    class vector_tests {
    private:
    public:
        vector_tests();

        static bool constructor_with_default_value_test();

        static bool addition_assignment_test();

        static bool subtraction_assignment_test();

        static bool multiplication_with_number_assignment_test();

        static bool division_by_number_assignment_test();

        static bool addition_test();

        static bool subtraction_test();

        static bool multiplication_with_number_test();

        static bool division_by_number_test();

    };
}

#endif //MAIN_CPP_VECTOR_TESTS_HPP
