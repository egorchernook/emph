#include <cassert>
#include <cmath>
#include <functional>
#include "one_step_methods_tests.hpp"

Numerical_methods::one_step_methods_tests::one_step_methods_tests() {
    std::cout << std::endl << "\ttesting one step ODE solving methods..." << std::endl;
    assert( Runge_Kutta_method_test() == true);
    std::cout << "\t\tRunge Kutta method works" << std::endl;
}

bool Numerical_methods::one_step_methods_tests::Runge_Kutta_method_test() {
    constexpr int precision_order = 4;
    vector<double> initial_state = { 0.0, 1.0 };
    Butcher_table<3> table{
            { { 1.0/6 , 0.0   , -1.0/6 },
                      { 1.0/12, 5.0/12,  0.0   },
                      { 0.5   , 1.0/3 ,  1.0/6 } },
            {   1.0/6 , 2.0/3 , 1.0/6    },
            {   0.0   , 1.0/2   , 1.0    }
    };
    implicit_Runge_Kutta_method< vector<double>, 3, precision_order > impl_rkm_method{
        initial_state,
        one_step_methods_tests::harmonic_oscillator_function };
    impl_rkm_method.set_Butcher_table( table);
    bool result = true;
    double time = 0.0;
    const double step = 0.000'1;
    for( int i =0; i < 100'000; ++i) {
        auto current_solution = impl_rkm_method.make_step(
                initial_state,
                step,
                time);
        const double x = current_solution.solution[0];
        const double y = current_solution.solution[1];
        if (std::abs(x - std::sin(time + step)) > step ||
            std::abs(y - std::cos(time + step)) > step ){
            result = false;
        }
        time += step;
        initial_state = { x, y};
    }
    return result;
}