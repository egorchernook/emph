#ifndef MAIN_CPP_BOUNDARY_VALUE_PROBLEM_HPP
#define MAIN_CPP_BOUNDARY_VALUE_PROBLEM_HPP

#include <functional>
#include <type_traits>
#include <cmath>

#include "ODE_solver.hpp"
#include "solution_t_concept.hpp"
#include "vector.hpp"

namespace Numerical_methods {

    template<
            solution_element_t solution_t,
            std::size_t precision_order = 4,
            template <solution_t_concept = vector<solution_t>, std::size_t = precision_order> class solver_t
                    = Adams_Moulton_method
                            >
    class shooting_method {
    public :
        using function_t = std::function<solution_t(double, const solution_t&, const solution_t&)>;
    private:
        solver_t<> ode_solver;
        solver_t<> ode_solver_for_s_derivative;
        const function_t function;
        const function_t function_derivative_y;
        const function_t function_derivative_u;
        solution_t current_angle;
        double step;
    public:
        solution_t get_current_angle(){
            return current_angle;
        }

        const struct boundary_conditions{
            const solution_t& left_value;
            const solution_t& right_value;
            const double left_border;
            const double right_border;
        } conditions;

        shooting_method( const function_t& function_,
                         const function_t& function_derivative_y_,
                         const function_t& function_derivative_u_,
                         const boundary_conditions& boundary_values,
                         const solution_t& initial_angle,
                         double step_size )
                         : function(function_),
                           function_derivative_y(function_derivative_y_),
                           function_derivative_u(function_derivative_u_),
                           conditions(boundary_values),
                           current_angle(initial_angle),
                           step(step_size) {

            ///only works for default template template parameter!!!
                //constexpr int precision_order = 4;
                Butcher_table<3> table{
                        {{1.0 / 6, 0.0, -1.0 / 6},
                                  {1.0 / 12, 5.0 / 12, 0.0},
                                           {0.5, 1.0 / 3, 1.0 / 6}},
                        {1.0 / 6, 2.0 / 3, 1.0 / 6},
                        {0.0,     1.0 / 2, 1.0}
                };
                implicit_Runge_Kutta_method<vector<double>, precision_order-1, precision_order> impl_rkm_method{};
                impl_rkm_method.set_Butcher_table(table);

                ode_solver.set_starting_method(impl_rkm_method);
                ode_solver.set_coefficients({3.0 / 8, 19.0 / 24, -5.0 / 24, 1.0 / 24});

                Butcher_table<3> another_table{
                        {{1.0 / 6, 0.0, -1.0 / 6},
                                  {1.0 / 12, 5.0 / 12, 0.0},
                                           {0.5, 1.0 / 3, 1.0 / 6}},
                        {1.0 / 6, 2.0 / 3, 1.0 / 6},
                        {0.0,     1.0 / 2, 1.0}
                };
                implicit_Runge_Kutta_method<vector<double>, precision_order-1, precision_order> another_impl_rkm_method{};
                impl_rkm_method.set_Butcher_table(table);

                ode_solver_for_s_derivative.set_starting_method(another_impl_rkm_method);
                ode_solver_for_s_derivative.set_coefficients({3.0 / 8, 19.0 / 24, -5.0 / 24, 1.0 / 24});
            ///only works for default template template parameter!!!

            ode_solver.set_initial_conditions( {conditions.left_value, current_angle} );
            ode_solver.set_function(
                    [this]( const vector<solution_t>& old_solution, double x ) -> vector<solution_t>
                    {
                        assert( old_solution.size() == 2);
                        return { old_solution[1], function(x, old_solution[0], old_solution[1]) };
                    } );
            ode_solver.set_time( conditions.left_border );

            if constexpr ( std::is_same_v<solution_t,double> ) {
                ode_solver_for_s_derivative.set_initial_conditions( {0.0,1.0});
            } else {
                solution_t first = conditions.left_value;
                for( auto &x : first) {
                    first = 0.0;
                }
                solution_t second = conditions.left_value;
                for( auto &x : second) {
                    second = 1.0;
                }
                ode_solver_for_s_derivative.set_initial_conditions( {first,second});
            }

            ode_solver_for_s_derivative.set_time( conditions.left_border);
        }

        return_t<vector<solution_t>> get_next_step() {
            assert( ode_solver.get_time() + step < conditions.right_border + std::numeric_limits<double>::epsilon()*10'000);
            auto result = ode_solver.get_next_step( step );
            if( ode_solver.buffer_is_filled() ) {
                auto last_elements = ode_solver.get_buffer_state();
                const auto current_solution = result.solution;
                ode_solver_for_s_derivative.set_function(
                        [this, &last_elements, &current_solution]
                                (const vector<solution_t> &old_solution, double x) -> vector<solution_t>
                                        {
                            assert(old_solution.size() == 2);
                            assert(std::abs(ode_solver.get_time() / step - x / step) <=
                                   std::numeric_limits<double>::epsilon() * 10'000 + ode_solver.size());

                            if (std::abs(x - ode_solver.get_time()) <
                                std::numeric_limits<double>::epsilon() * 10'000) {
                                return {old_solution[1],
                                        function_derivative_y(x, current_solution[0], current_solution[1]) *
                                        old_solution[0] +
                                        function_derivative_u(x, current_solution[0], current_solution[1]) *
                                        old_solution[1]};
                            } else {
                                const int idx = std::round(std::abs(ode_solver.get_time() / step - x / step)) - 1;
                                return {old_solution[1],
                                        function_derivative_y(x, last_elements[idx][0], last_elements[idx][1]) *
                                        old_solution[0] +
                                        function_derivative_u(x, last_elements[idx][0], last_elements[idx][1]) *
                                        old_solution[1]};
                            }
                        }
                );
            } else {
                auto last_element = ode_solver.get_last_element();
                ode_solver_for_s_derivative.set_function(
                        [this, &last_element]
                                (const vector<solution_t> &old_solution, double x) -> vector<solution_t>
                        {
                            assert(old_solution.size() == 2);
                            assert(std::abs(ode_solver.get_time() / step - x / step) <=
                                   std::numeric_limits<double>::epsilon() * 10'000 + ode_solver.size());

                            return {    old_solution[1],
                                        function_derivative_y( x, last_element[0], last_element[1]) *
                                        old_solution[0] +
                                        function_derivative_u( x, last_element[0], last_element[1]) *
                                        old_solution[1]};
                        }
                );
            }
            auto need_to_make_a_step = ode_solver_for_s_derivative.get_next_step( step);
            return result;
        }

        solution_t get_new_angle(){
            return current_angle
                    - inverse(ode_solver_for_s_derivative.get_buffer_state()[0])
                    * ode_solver.get_buffer_state()[0]
                    - conditions.right_value;
        }

        double get_solution_error(){
            return ode_solver.get_buffer_state()[0] - conditions.right_value;
        }

        void restart( const solution_t& angle ){
            current_angle = angle;
            auto new_ode_solver = decltype(ode_solver)();
            auto new_ode_s_solver = decltype(ode_solver_for_s_derivative)();
            ode_solver = new_ode_solver;
            ode_solver_for_s_derivative = new_ode_s_solver;

            ode_solver.set_initial_conditions( {conditions.left_value, current_angle} );
            ode_solver.set_function(
                    [this]( const vector<solution_t>& old_solution, double x ) -> vector<solution_t>
                    {
                        assert( old_solution.size() == 2);
                        return { old_solution[1], function(x, old_solution[0], old_solution[1]) };
                    } );
            ode_solver.set_time( conditions.left_border );

            if constexpr ( std::is_same_v<solution_t,double> ) {
                ode_solver_for_s_derivative.set_initial_conditions( {0.0,1.0});
            } else {
                solution_t first = conditions.left_value;
                for( auto &x : first) {
                    first = 0.0;
                }
                solution_t second = conditions.left_value;
                for( auto &x : second) {
                    second = 1.0;
                }
                ode_solver_for_s_derivative.set_initial_conditions( {first,second});
            }

            ode_solver_for_s_derivative.set_time( conditions.left_border);
        }
    };
}
#endif //MAIN_CPP_BOUNDARY_VALUE_PROBLEM_HPP