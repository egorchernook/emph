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
#include "boundary_value_problem.hpp"
#include "heat_equation_2d_solver.hpp"

#include "tests.hpp"

void harmonic(){
    constexpr int steps_amount = 10'000;
    constexpr double period = 10.0 * 3.141'592'654;
    constexpr double step = period / steps_amount;
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
            harmonic_oscillator_function};
    impl_rkm_method.set_Butcher_table(table);
    Adams_Moulton_method<vector<double>, precision_order> impl_msm_solver{
            initial_state,
            harmonic_oscillator_function
    };

    impl_msm_solver.set_starting_method(impl_rkm_method);
    impl_msm_solver.set_coefficients({3.0 / 8, 19.0 / 24, -5.0 / 24, 1.0 / 24});

    std::ofstream output("../../results/test/test1.dat", std::ios_base::trunc);
    for( int i = 0; i < steps_amount; ++i ) {
        const auto result = impl_msm_solver.get_next_step( step );

        output << i*step << "\t" << result.solution[0] << "\t" << result.solution[1] << "\t\t" << std::sin( (i + 1)*step ) << "\t" << std::cos( (i + 1)*step) << "\n";
    }

    output.flush();
    output.close();
}
void first_task() {
        constexpr int steps_amount = 100'000;
        constexpr double period = 100.0 * 3.141'592'654;
        constexpr double step = period / steps_amount;
        std::string filename = "../../results/first/res_";
        using namespace Numerical_methods;

        const int precision_order = 4;
        Butcher_table<3> table{
            {{1.0 / 6, 0.0, -1.0 / 6},
                {1.0 / 12, 5.0 / 12, 0.0},
                {0.5, 1.0 / 3, 1.0 / 6}},
            {1.0 / 6, 2.0 / 3, 1.0 / 6},
            {0.0,     1.0 / 2, 1.0}
        };
        std::vector<double> coefs{3.0 / 8, 19.0 / 24, -5.0 / 24, 1.0 / 24};
        constexpr double omega2 = 1.0;
        constexpr int omegas_amount = 5;
        /*
        for(int o1 = 0; o1 < omegas_amount; ++o1 ){
            const double omega1 = omega2 / std::pow( 10.0, 2 * o1 + 1);
            std::cout << "o1 : " << o1 << std::endl;
            for( int o3 = 0; o3 < omegas_amount; ++o3 ) {
                const double omega3 = omega2 / std::pow(10.0, 2 * o3 + 1);
                std::cout << "o3 : " << o3 << std::endl;
                */
                const double omega1 = 0.0001;
                const double omega3 = 0.0001;
                vector<double> initial_state = {omega1, omega2, omega3, 0.0, 0.0, 0.0};
                const std::string filename_omega = std::to_string(omega1) + "_" +
                                                   std::to_string(omega2) + "_" +
                                                   std::to_string(omega3) + "_";
                /*
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
                */
                const double i2 = 0.5;
                const double i3 = 0.1;
                const double energy = 1.0 * omega1 * omega1 + i2 * omega2 * omega2 + i3 * omega1 * omega1;
                auto function = [&i2, &i3](const vector<double> &old_solution, double time) noexcept
                        -> vector<double> {
                    const auto first = (i2 - i3) * old_solution[1] * old_solution[2];
                    const auto second = (i3 - 1.0) / i2 * old_solution[0] * old_solution[2];
                    const auto third = (1.0 - i2) / i3 * old_solution[0] * old_solution[1];
                    return {first, second, third, old_solution[0], old_solution[1], old_solution[2]};
                };
                implicit_Runge_Kutta_method<vector<double>, 3, precision_order> impl_rkm_method{
                        initial_state,
                        function};
                impl_rkm_method.set_Butcher_table(table);
                Adams_Moulton_method<vector<double>, precision_order> impl_msm_solver{
                        initial_state,
                        function
                };
                impl_msm_solver.set_starting_method(impl_rkm_method);
                impl_msm_solver.set_coefficients(coefs);

                std::ofstream output(
                        filename + filename_omega + /*filename_i2*/ "0.5" + "_" + /*filename_i3*/ "0.135" + ".dat",
                        std::ios_base::trunc);

                for (int i = 0; i < steps_amount; ++i) {
                    const auto result = impl_msm_solver.get_next_step(step);

                    output << i * step
                           << "\t" << result.solution[0] / omega1
                           << "\t" << result.solution[1] / omega2
                           << "\t" << result.solution[2] / omega3
                           << "\t"
                           << "\t" << result.solution[3]
                           << "\t" << result.solution[4]
                           << "\t" << result.solution[5]
                           << std::endl;
                }

                output.flush();
                output.close();
                /*}
            }    */
            //}
        //}
};
void second_task_test1(){
    using namespace Numerical_methods;

    constexpr int precision_order = 4;
    Butcher_table<precision_order-1> table{
            {{1.0 / 6, 0.0, -1.0 / 6},
                      {1.0 / 12, 5.0 / 12, 0.0},
                               {0.5, 1.0 / 3, 1.0 / 6}},
            {1.0 / 6, 2.0 / 3, 1.0 / 6},
            {0.0,     1.0 / 2, 1.0}
    };
    std::vector<double> coefficients = {3.0 / 8, 19.0 / 24, -5.0 / 24, 1.0 / 24};
    implicit_Runge_Kutta_method<vector<double>, precision_order-1, precision_order> impl_rkm_method{};
    impl_rkm_method.set_Butcher_table(table);

    Adams_Moulton_method< vector<double>, precision_order> impl_msm_solver{};
    impl_msm_solver.set_starting_method( impl_rkm_method );
    impl_msm_solver.set_coefficients( coefficients );

    implicit_Runge_Kutta_method<vector<double>, precision_order-1, precision_order> another_impl_rkm_method{};
    another_impl_rkm_method.set_Butcher_table(table);

    Adams_Moulton_method< vector<double>, precision_order> another_impl_msm_solver{};
    another_impl_msm_solver.set_starting_method(another_impl_rkm_method);
    another_impl_msm_solver.set_coefficients( coefficients );

    constexpr double c1 = 1.0;
    constexpr double c2 = 1.0;
    const auto exact_result = [&c1,&c2](double x) -> double {
        return c1*sin(x) + c2*cos(x);
    };
    const auto exact_result_derivative = [&c1,&c2](double x) -> double {
        return c1*cos(x) - c2*sin(x);
    };

    const boundary_conditions_for_shooting_method<double> conditions{1.0, exact_result(10), 0.0, 10.0};

    constexpr double step = 0.000'1;
    shooting_method<double> solver(
            [](double x, const double& y, const double& y1) -> double {
                return -y;
            },
            [](double x, const double& y, const double& y1) -> double {
                return -1.0;
            },
            [](double x, const double& y, const double& y1) -> double {
                return 0.0;
            },
            conditions,
            3.0,
            step,
            impl_msm_solver,
            another_impl_msm_solver
    );

    const int steps_amount = (conditions.right_border - conditions.left_border)/step;
    int counter = 1;

    std::vector<std::vector<double>> data{};
    do {
        data.clear();
        for (int i = 0; i <= steps_amount; ++i) {
            double y{0.0};
            double y1{0.0};
            double Y{0.0};
            double U{0.0};
            if( i == 0 ) {
                y = solver.conditions.left_value;
                y1 = solver.get_current_angle();
                Y = 0.0;
                U = 1.0;
            } else {
                const auto result = solver.get_next_step();
                y = result[0].solution[0];
                y1 = result[0].solution[1];
                Y = result[1].solution[0];
                U = result[1].solution[1];
            }
            data.push_back( {
                            i * step,
                            y,
                            y1,
                            Y,
                            U,
                            exact_result(i * step),
                            exact_result_derivative( i * step)
            });
        }
        if( solver.is_accuracy_reached( step*step*step ) ){
            break;
        } else {
            std::cout << solver.get_solution_error() << "\t" << solver.get_current_angle() << std::endl;
            solver.restart( solver.get_new_angle() );
            counter++;
        }
    } while ( true );
    std::ofstream output("../../results/test/test_second1_test.dat", std::ios_base::trunc);
    for( auto& x : data) {
        output << x[0]
               << "\t" << x[1] << "\t" << x[2] << "\t" << x[3] << "\t" << x[4]
               << "\t\t" << x[5]
               << "\t" << x[6]
               << "\n";
    }
    output.flush();
    output.close();
    std::cout << "shots : " << counter << "\t" << "angle : " << solver.get_current_angle() << std::endl;
};
void second_task_test2(){
    using namespace Numerical_methods;

    constexpr int precision_order = 4;
    Butcher_table<precision_order-1> table{
            {{1.0 / 6, 0.0, -1.0 / 6},
                      {1.0 / 12, 5.0 / 12, 0.0},
                               {0.5, 1.0 / 3, 1.0 / 6}},
            {1.0 / 6, 2.0 / 3, 1.0 / 6},
            {0.0,     1.0 / 2, 1.0}
    };
    std::vector<double> coefficients = {3.0 / 8, 19.0 / 24, -5.0 / 24, 1.0 / 24};
    implicit_Runge_Kutta_method<vector<double>, precision_order-1, precision_order> impl_rkm_method{};
    impl_rkm_method.set_Butcher_table(table);

    Adams_Moulton_method< vector<double>, precision_order> impl_msm_solver{};
    impl_msm_solver.set_starting_method( impl_rkm_method );
    impl_msm_solver.set_coefficients( coefficients );

    implicit_Runge_Kutta_method<vector<double>, precision_order-1, precision_order> another_impl_rkm_method{};
    another_impl_rkm_method.set_Butcher_table(table);

    Adams_Moulton_method< vector<double>, precision_order> another_impl_msm_solver{};
    another_impl_msm_solver.set_starting_method(another_impl_rkm_method);
    another_impl_msm_solver.set_coefficients( coefficients );

    constexpr double c1 = 1.0;
    constexpr double c2 = 1.0;
    auto ch = [](double x) -> double {
        return 0.5 * (std::exp(x) + std::exp(-x));
    };
    auto sh = [](double x) -> double {
        return 0.5 * (std::exp(x) - std::exp(-x));
    };
    const auto exact_result = [ch,sh](double x) -> double {
        return c1*sh(x) + c2*ch(x);
    };
    const auto exact_result_derivative = [&c1,&c2,ch,sh](double x) -> double {
        return c1*ch(x) + c2*sh(x);
    };

    const boundary_conditions_for_shooting_method<double> conditions{1.0, exact_result(10), 0.0, 10.0};

    constexpr double step = 0.000'1;
    shooting_method<double> solver(
            [](double x, const double& y, const double& y1) -> double {
                return y;
            },
            [](double x, const double& y, const double& y1) -> double {
                return 1.0;
            },
            [](double x, const double& y, const double& y1) -> double {
                return 0.0;
            },
            conditions,
            3.0,
            step,
            impl_msm_solver,
            another_impl_msm_solver
    );

    const int steps_amount = (conditions.right_border - conditions.left_border)/step;
    int counter = 1;

    std::vector<std::vector<double>> data{};
    do {
        data.clear();
        for (int i = 0; i <= steps_amount; ++i) {
            double y{0.0};
            double y1{0.0};
            double Y{0.0};
            double U{0.0};
            if( i == 0 ) {
                y = solver.conditions.left_value;
                y1 = solver.get_current_angle();
                Y = 0.0;
                U = 1.0;
            } else {
                const auto result = solver.get_next_step();
                y = result[0].solution[0];
                y1 = result[0].solution[1];
                Y = result[1].solution[0];
                U = result[1].solution[1];
            }
            data.push_back( {
                                    i * step,
                                    y,
                                    y1,
                                    Y,
                                    U,
                                    exact_result(i * step),
                                    exact_result_derivative( i * step)
                            });
        }
        if( solver.is_accuracy_reached( step*step ) ){
            break;
        } else {
            std::cout << solver.get_solution_error() << "\t" << solver.get_current_angle() << std::endl;
            solver.restart( solver.get_new_angle() );
            counter++;
        }
    } while ( true );
    std::ofstream output("../../results/test/test_second2.dat", std::ios_base::trunc);
    for( auto& x : data) {
        output << x[0]
               << "\t" << x[1] << "\t" << x[2] << "\t" << x[3] << "\t" << x[4]
               << "\t\t" << x[5]
               << "\t" << x[6]
               << "\n";
    }
    output.flush();
    output.close();
    std::cout << "shots : " << counter << "\t" << "angle : " << solver.get_current_angle() << std::endl;
};
void second_task() {
    using namespace Numerical_methods;

    constexpr int precision_order = 4;
    constexpr int table_size = 3;
    Butcher_table<table_size> table{
            {{1.0 / 6, 0.0, -1.0 / 6},
                      {1.0 / 12, 5.0 / 12, 0.0},
                               {0.5, 1.0 / 3, 1.0 / 6}},
            {1.0 / 6, 2.0 / 3, 1.0 / 6},
            {0.0,     1.0 / 2, 1.0}
    };
    std::vector<double> coefficients = {3.0 / 8, 19.0 / 24, -5.0 / 24, 1.0 / 24};
    implicit_Runge_Kutta_method<vector<double>, table_size, precision_order> impl_rkm_method{};
    impl_rkm_method.set_Butcher_table(table);

    Adams_Moulton_method<vector<double>, precision_order> impl_msm_solver{};
    impl_msm_solver.set_starting_method(impl_rkm_method);
    impl_msm_solver.set_coefficients(coefficients);

    implicit_Runge_Kutta_method<vector<double>, table_size, precision_order> another_impl_rkm_method{};
    another_impl_rkm_method.set_Butcher_table(table);

    Adams_Moulton_method<vector<double>, precision_order> another_impl_msm_solver{};
    another_impl_msm_solver.set_starting_method(another_impl_rkm_method);
    another_impl_msm_solver.set_coefficients(coefficients);

    const auto lambda = [](double x, double y) -> double {
        return y*y*y;
    };
    const auto f = [](double x, double y) -> double {
        return -x*y*y*y*y;
    };
    const auto q = [lambda](double x, double y, double y1) -> double {
        return -lambda(x, y) * y1;
    };

    const boundary_conditions_for_shooting_method<double> conditions{0.0, 1.0, 0.0, 1.0};

    constexpr double step = 0.001;
    constexpr double epsilon = step / 100;
    shooting_method<double> solver(
            [epsilon](double x, const double &y, const double &y1) -> double {
                double result{0.0};
                if ( std::abs(y) < epsilon ){
                    result = x * y + 3.0 * y1 * y1 / epsilon;
                } else {
                    result = x * y + 3.0 * y1 * y1 / y;
                }
                return result;
            },
            [epsilon](double x, const double &y, const double &y1) -> double {
                double result{0.0};
                if ( std::abs(y * y) < epsilon ){
                    result = x + 3.0 * y1 * y1 / epsilon;
                } else {
                    result = x * y + 3.0 * y1 * y1 / (y * y);
                }
                return result;
            },
            [epsilon](double x, const double &y, const double &y1) -> double {
                double result{0.0};
                if ( std::abs(y) < epsilon ) {
                    result = -6.0 * y1 / epsilon;
                } else {
                    result = -6.0 * y1 / y;
                }
                return result;
            },
            conditions,
            3.0,
            step,
            impl_msm_solver,
            another_impl_msm_solver
    );
    const int steps_amount = (conditions.right_border - conditions.left_border) / step;
    int counter = 1;

    std::ofstream shooting_parameters("../../results/second/shooting_parameters.dat", std::ios_base::trunc);
    do {
        std::ofstream output("../../results/second/second" + std::to_string(counter) + ".dat", std::ios_base::trunc);
        double y{0.0};
        double y1{0.0};
        double Y{0.0};
        double U{0.0};
        for (int i = 0; i <= steps_amount; ++i) {
            const double x = i * step;
            if (i == 0) {
                y = solver.conditions.left_value;
                y1 = solver.get_current_angle();
                Y = 0.0;
                U = 1.0;
            } else {
                const auto result = solver.get_next_step();
                y = result[0].solution[0];
                y1 = result[0].solution[1];
                Y = result[1].solution[0];
                U = result[1].solution[1];
            }
            output  << x << "\t" << y << "\t" << y1 << "\t"
                    << q( x, y, y1) << "\t" << lambda( x, y) << "\t" << f( x, y) << "\t\t"
                    << Y << "\t" << U
                    << std::endl;
        }

        output.flush();
        output.close();

        if (solver.is_accuracy_reached(step/100)) {
            shooting_parameters << counter << "\t" << solver.get_current_angle() << std::endl;
            break;
        } else {
            shooting_parameters << counter << "\t" << solver.get_current_angle() << std::endl;
            std::cout << solver.get_solution_error() << "\t" << solver.get_current_angle() << std::endl;
            solver.restart(solver.get_new_angle());
            counter++;
        }
    } while (false);

    shooting_parameters.flush();
    shooting_parameters.close();

    std::cout << "shots : " << counter << "\t" << "final angle : " << solver.get_current_angle() << std::endl;
};

void third_test(){
    constexpr int steps_amount = 10'000;
    constexpr double step = 0.01;

    using namespace Numerical_methods;

    constexpr int precision_order = 4;
    constexpr int table_size = 3;
    Butcher_table<table_size> table{
            {{1.0 / 6, 0.0, -1.0 / 6},
                      {1.0 / 12, 5.0 / 12, 0.0},
                               {0.5, 1.0 / 3, 1.0 / 6}},
            {1.0 / 6, 2.0 / 3, 1.0 / 6},
            {0.0,     1.0 / 2, 1.0}
    };
    std::vector<double> coefficients = {3.0 / 8, 19.0 / 24, -5.0 / 24, 1.0 / 24};
    implicit_Runge_Kutta_method<matrix<double>, table_size, precision_order> impl_rkm_method{};
    impl_rkm_method.set_Butcher_table(table);

    Adams_Moulton_method<matrix<double>, precision_order> impl_msm_solver{};
    impl_msm_solver.set_starting_method(impl_rkm_method);
    impl_msm_solver.set_coefficients(coefficients);

    const first_and_second_conditions<double> conditions{
            {[](double, double) -> double {
                return 0;
            },
                    [](double, double) -> double {
                        return 0;
                    }
            },
            {[](double, double) -> double {
                return 0;
            },
                    [](double, double) -> double {
                        return 0;
                    }
            },
            {0.0, 0.0},
            {1.0, 1.0}
    };

    heat_equation_2d_solver<precision_order> solver{
            [](const vector<double> &r) -> double {
                return std::exp(r[0]);
            },
            [](const vector<double>&) -> double {
                return 0.25;
            },
            [](const vector<double>&) -> double {
                return 1.0;
            },
            [](const vector<double> &, double) -> double {
                return 0.0;
            },
            &conditions,
            {0.05, 0.05},
            impl_msm_solver,
            0.0,
            step
    };


    std::ofstream output("../../results/test/test3/test3.dat", std::ios_base::trunc);
    for( int i = 0; i < steps_amount; ++i ) {
        std::cout << "Time : " << i * step << std::endl;
        const auto result = solver.get_next_step();
        const auto layer = result.solution;
        for( int j = 0; j < layer.height(); ++j ) {
            for( int k = 0; k < layer.width(); ++k ) {
                output  << i * step << "\t"
                        << j * solver.lattice_steps_sizes[0] << "\t"
                        << k * solver.lattice_steps_sizes[1] << "\t"
                        << layer[j][k] << std::endl;
            }
        }
    }

    output.flush();
    output.close();
}

int main() {
    //Numerical_methods::test();
    //harmonic();
    //first_task();
    //second_task_test1();
    //second_task_test2();
    //second_task();
    third_test();
    return 0;
};