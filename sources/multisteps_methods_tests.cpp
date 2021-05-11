#include <cassert>
#include <cmath>
#include "multisteps_methods_tests.hpp"
#include "one_step_methods.hpp"
#include "ring_buffer.hpp"

Numerical_methods::multisteps_methods_tests::multisteps_methods_tests() {
    std::cout << std::endl << "\ttesting multisteps ODE solving methods..." << std::endl;
    assert( Adams_Moulton_method_test() );
    std::cout << "\t\tAdams Moulton method works" << std::endl;
}

bool Numerical_methods::multisteps_methods_tests::Adams_Moulton_method_test() {
    vector<double> initial_state = { 0.0, 1.0 };
    constexpr int precision_order = 4;
    Butcher_table<3> table{
            { { 1.0/6 , 0.0   , -1.0/6 },
                        { 1.0/12, 5.0/12,  0.0   },
                                  { 0.5   , 1.0/3 ,  1.0/6 } },
            {   1.0/6 , 2.0/3 , 1.0/6    },
            {   0.0   , 1.0/2   , 1.0    }
    };
    implicit_Runge_Kutta_method< vector<double>, 3, precision_order > impl_rkm_method{
        initial_state,
        multisteps_methods_tests::harmonic_oscillator_function,
        };
    impl_rkm_method.set_Butcher_table(table);
    Adams_Moulton_method< vector<double>, precision_order> impl_msm_solver{
        initial_state,
        multisteps_methods_tests::harmonic_oscillator_function
        };
    impl_msm_solver.set_starting_method( impl_rkm_method );
    impl_msm_solver.set_coefficients( {3.0/8, 19.0/24, -5.0/24, 1.0/24} );
    ///
    implicit_Runge_Kutta_method< vector<double>, 3, precision_order > another_impl_rkm_method{};
    another_impl_rkm_method.set_Butcher_table(table);
    another_impl_rkm_method.set_function( multisteps_methods_tests::harmonic_oscillator_function);
    auto initial_state1 = initial_state;
    another_impl_rkm_method.set_initial_conditions( initial_state1);
    another_impl_rkm_method.set_time( 0.0);

    impl_rkm_method.set_Butcher_table(table);
    Adams_Moulton_method< vector<double>, precision_order> another_impl_msm_solver{};
    ///работает только в таком порядке!!!
    another_impl_msm_solver.set_initial_conditions(initial_state1);
    another_impl_msm_solver.set_function(multisteps_methods_tests::harmonic_oscillator_function);
    another_impl_msm_solver.set_starting_method( impl_rkm_method );
    another_impl_msm_solver.set_time(0.0);
    another_impl_msm_solver.set_coefficients( {3.0/8, 19.0/24, -5.0/24, 1.0/24} );

    bool result = true;
    double time = 0.0;
    const double step = 0.000'1;
    ring_buffer< vector<double>, precision_order - 1 > buffer;
    ring_buffer< vector<double>, precision_order - 1 > buffer1;
    for(int i = 0; i < 3; ++i){
        auto current_solution = impl_rkm_method.make_step(
                initial_state,
                step,
                time + step );
        auto another_current_solution = another_impl_rkm_method.make_step(
                initial_state,
                step,
                time + step );
        buffer.push( current_solution.solution );
        buffer1.push( another_current_solution.solution );
        time += step;
        initial_state = current_solution.solution;
        initial_state1 = another_current_solution.solution;
    }
    for(int i = 0; i < 100'000; ++i) {
        time += step;
        const auto temp = buffer.current_state();
        const auto temp1 = buffer1.current_state();
        const std::vector< vector<double>> last_steps{ temp.rbegin(), temp.rend()};
        const std::vector< vector<double>> last_steps1{ temp1.rbegin(), temp1.rend()};
        auto msm_solution = impl_msm_solver.make_step(
                last_steps,
                step,
                time);
        auto msm_solution1 = another_impl_msm_solver.make_step(
                last_steps1,
                step,
                time);
        const double x = msm_solution.solution[0];
        const double y = msm_solution.solution[1];
        const double x1 = msm_solution1.solution[0];
        const double y1 = msm_solution1.solution[1];
        if (std::abs(x - std::sin( time)) > step ||
            std::abs(y - std::cos( time)) > step ||
            std::abs(x1 - std::sin( time)) > step ||
            std::abs(y1 - std::cos( time)) > step) {
            result = false;
        }
        buffer.push( { x, y} );
        buffer1.push( { x1, y1} );
    }
    return result;
}