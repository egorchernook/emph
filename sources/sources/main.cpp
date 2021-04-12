#include <iostream>
#include <fstream>
#include <algorithm>

#include "solution_t_concept.hpp"
#include "ring_buffer.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "one_step_methods.hpp"
#include "multisteps_methods.hpp"

#include "tests.hpp"

constexpr int steps_amount = 100'000;
constexpr double step = 2*3.1415 / steps_amount;
using namespace Numerical_methods;
int main(){

    auto harmonic_oscillator_function = [](const vector<double>& old_solution, double time) noexcept
            -> vector<double> {
        return { old_solution[1], -old_solution[0] };
    };

    vector<double> initial_state = {0.0, 1.0};
    constexpr int precision_order = 4;
    Butcher_table<3> table{
        {{1.0 / 6, 0.0, -1.0 / 6},
                 {1.0 / 12, 5.0 / 12, 0.0},
                 {0.5, 1.0 / 3, 1.0 / 6}},
         {1.0 / 6, 2.0 / 3, 1.0 / 6},
         {0.0,     1.0 / 2, 1.0}
    };
    implicit_Runge_Kutta_method<vector<double>, 3, precision_order> impl_rkm_method{
        initial_state,
        harmonic_oscillator_function,};
    impl_rkm_method.set_Butcher_table(table);
    Adams_Moulton_method<vector<double>, precision_order> impl_msm_solver{
        initial_state,
        harmonic_oscillator_function
    };

    impl_msm_solver.set_starting_method(impl_rkm_method);
    impl_msm_solver.set_coefficients({3.0 / 8, 19.0 / 24, -5.0 / 24, 1.0 / 24});

    std::ofstream output("../../results/test.dat", std::ios_base::trunc);
    for( int i = 0; i < steps_amount; ++i ) {
        const auto result = impl_msm_solver.get_next_step( step );

        output << i*step << "\t" << result.solution[0] << "\t" << result.solution[1] << "\t\t" << std::sin( i*step ) << "\t" << std::cos( i*step) << "\n";
    }

    output.flush();
    output.close();
    return 0;
}