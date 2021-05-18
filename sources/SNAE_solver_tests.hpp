#ifndef MAIN_CPP_SNAE_SOLVER_TESTS_HPP
#define MAIN_CPP_SNAE_SOLVER_TESTS_HPP

#include <cmath>
#include "vector.hpp"
#include "SNAE_solver.hpp"

namespace Numerical_methods {
    class SNAE_solver_tests {
    private:
        /* system : x + y = 1; 5y - x = 2; solution : x = 1/2; y = 1/2; */
        static vector<double> test_equation(const vector<double> &input_data) noexcept {
            assert( input_data.size() == 2);
            return { 1 - input_data[1], 2.0/5 + input_data[0] / 5 };
        }
        static vector<double> test_equation_derivative_1(const vector<double> &input_data) noexcept {
            assert( input_data.size() == 2);
            return { 1.0, 1.0};
        }
        static vector<double> test_equation_derivative_2(const vector<double> &input_data) noexcept {
            assert( input_data.size() == 2);
            return { -1.0/5.0, 1.0};
        }
    public:
        SNAE_solver_tests();
        static bool test_fixed_point_iterations_method();
        static bool test_Seidels_method();
        static bool test_Newtons_method();
    };
}

#endif //MAIN_CPP_SNAE_SOLVER_TESTS_HPP
