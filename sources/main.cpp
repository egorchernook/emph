#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>

#include "solution_t_concept.hpp"
#include "ring_buffer.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "one_step_methods.hpp"
#include "multisteps_methods.hpp"

#include "tests.hpp"

void harmonic(){
    constexpr int steps_amount = 100'000;
    constexpr double step = 50*3.1415 / steps_amount;
    using namespace Numerical_methods;
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
}

int main(){

    constexpr int steps_amount = 100'000;
    constexpr double period = 100.0 * 3.141'592'654;
    constexpr double step = period / steps_amount;
    std::string filename = "../../results/res_";
    using namespace Numerical_methods;

    constexpr int precision_order = 4;
    Butcher_table<3> table{
            {{1.0 / 6, 0.0, -1.0 / 6},
             {1.0 / 12, 5.0 / 12, 0.0},
             {0.5, 1.0 / 3, 1.0 / 6}},
            {1.0 / 6, 2.0 / 3, 1.0 / 6},
            {0.0,     1.0 / 2, 1.0}
    };
    std::vector<double> coefs{3.0 / 8, 19.0 / 24, -5.0 / 24, 1.0 / 24};
    constexpr double omega2 = 1.0;
        const double omega1 = 0.0001;
        vector<double> initial_state = { omega1, omega2, omega1,  0.0, 0.0, 0.0 };
        const std::string filename_omega = std::to_string(omega1) + "_" +
                                           std::to_string(omega2) + "_";
        const int i_steps = 10;
        for (int i2_step = 1; i2_step < i_steps; ++i2_step) {
            const double i2 = 1.0 * i2_step / i_steps;
            std::cout << "i2 : " << i2 << std::endl;
            const std::string filename_i2 = std::to_string(i2_step);
            const double i3_step_size = i2 / (i_steps + 1.0);
            for (int i3_step = 0; i3_step < i_steps; ++i3_step) {
                const double i3 = i3_step_size * (i3_step + 1);
                std::cout << "i3 : " << i3 << std::endl;
                const std::string filename_i3 = std::to_string(i3_step);
                auto function = [&i2,&i3](const vector<double> &old_solution, double time) noexcept
                        -> vector<double> {
                    const auto first = (i2 - i3) * old_solution[1] * old_solution[2];
                    const auto second = (i3 - 1.0) / i2 * old_solution[0] * old_solution[2];
                    const auto third = (1.0 - i2) / i3 * old_solution[0] * old_solution[1];
                    return { first, second, third, old_solution[0],  old_solution[1],  old_solution[2]};
                };
                implicit_Runge_Kutta_method<vector<double>, 3, precision_order> impl_rkm_method{
                        initial_state,
                        function,};
                impl_rkm_method.set_Butcher_table(table);
                Adams_Moulton_method<vector<double>, precision_order> impl_msm_solver{
                        initial_state,
                        function
                };
                impl_msm_solver.set_starting_method(impl_rkm_method);
                impl_msm_solver.set_coefficients(coefs);

                std::ofstream output(filename + filename_omega + filename_i2 + "_" + filename_i3 + ".dat", std::ios_base::trunc);

                for (int i = 0; i < steps_amount; ++i) {
                    const auto result = impl_msm_solver.get_next_step(step);

                    output << i * step
                           << "\t" << result.solution[0] /*/ omega1*/
                           << "\t" << result.solution[1]
                           << "\t" << result.solution[2] /*/ omega1*/
                           << "\t"
                           << "\t" << result.solution[3]
                           << "\t" << result.solution[4]
                           << "\t" << result.solution[5]
                           << std::endl;
                }

                output.flush();
                output.close();
            }
        }
    //}
    return 0;
}