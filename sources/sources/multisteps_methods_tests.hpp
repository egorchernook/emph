#ifndef MAIN_CPP_MULTISTEPS_METHODS_TESTS_HPP
#define MAIN_CPP_MULTISTEPS_METHODS_TESTS_HPP

#include "matrix.hpp"
#include "vector.hpp"
#include "multisteps_methods.hpp"
namespace Numerical_methods {
    class multisteps_methods_tests {
    private:
        static constexpr double omega_squared = 1.0;
        static vector<double> harmonic_oscillator_function( const vector<double>& old_solution, double time ) {
            assert( old_solution.size() == 2 );
            return { old_solution[1], -omega_squared * old_solution[0] };
        }
    public:
        multisteps_methods_tests();

        static bool Adams_Moulton_method_test();
    };
}

#endif //MAIN_CPP_MULTISTEPS_METHODS_TESTS_HPP
