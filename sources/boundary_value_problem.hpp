#ifndef MAIN_CPP_BOUNDARY_VALUE_PROBLEM_HPP
#define MAIN_CPP_BOUNDARY_VALUE_PROBLEM_HPP

#include <functional>
#include <type_traits>
#include <cmath>

#include "ODE_solver.hpp"
#include "solution_t_concept.hpp"
#include "vector.hpp"

namespace Numerical_methods {

    template<solution_element_t solution_t>
    struct boundary_conditions{
        const solution_t left_value;
        const solution_t right_value;
        const double left_border;
        const double right_border;
    };

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
        const boundary_conditions<solution_t>& conditions;

        shooting_method( const function_t& function_,
                         const function_t& function_derivative_y_,
                         const function_t& function_derivative_u_,
                         const boundary_conditions<solution_t>& boundary_values,
                         const solution_t& initial_angle,
                         double step_size,
                         solver_t<>& main_solver,
                         solver_t<>& additional_solver )
                         : function(function_),
                           function_derivative_y(function_derivative_y_),
                           function_derivative_u(function_derivative_u_),
                           conditions(boundary_values),
                           current_angle(initial_angle),
                           step(step_size),
                           ode_solver(main_solver),
                           ode_solver_for_s_derivative(additional_solver){


            ode_solver.set_initial_conditions( {conditions.left_value, current_angle} );
            auto func = [this]( const vector<solution_t>& old_solution, double x ) -> vector<solution_t>
            {
                assert( old_solution.size() == 2);
                return { old_solution[1], this->function(x, old_solution[0], old_solution[1]) };
            };
            ode_solver.set_function( func );
            ode_solver.set_time( conditions.left_border );
            if constexpr ( precision_order != 1){
                ode_solver.starting_method->set_initial_conditions( {conditions.left_value, current_angle} );
                ode_solver.starting_method->set_function( func );
                ode_solver.starting_method->set_time( conditions.left_border );
            }

            if constexpr ( std::is_same_v<solution_t,double> ) {
                ode_solver_for_s_derivative.set_initial_conditions( {0.0, 1.0});
                if constexpr ( precision_order != 1) {
                    ode_solver_for_s_derivative.starting_method->set_initial_conditions({0.0, 1.0});
                }
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
                if constexpr ( precision_order != 1){
                    ode_solver_for_s_derivative.starting_method->set_initial_conditions( {first,second} );
                }
            }
            ode_solver_for_s_derivative.set_time( conditions.left_border);
            if constexpr ( precision_order != 1){
                ode_solver_for_s_derivative.starting_method->set_time( conditions.left_border);
            }
        }

        std::vector<return_t<vector<solution_t>>> get_next_step() {

            bool have_starting_method = false;
            if constexpr( precision_order != 1){
                have_starting_method = true;
            }

            const auto time = ode_solver.get_time();
            assert( ode_solver.get_time() + step < conditions.right_border + std::numeric_limits<double>::epsilon()*10'000);

            if( ode_solver.buffer_is_filled() ) {
                const auto last_elements = ode_solver.get_buffer_state();

                ode_solver_for_s_derivative.set_function
                (
                        [this, last_elements, time]
                                (const vector<solution_t> &old_solution, double x) -> vector<solution_t>
                                        {
                            assert(old_solution.size() == 2);
                            assert(std::abs( time - x) / step <
                                   std::numeric_limits<double>::epsilon() * 10'000 + ode_solver.size());
                            const int idx = std::round(std::abs(time / step - x / step));
                            return {    old_solution[1],
                                        function_derivative_y(x, last_elements[idx][0], last_elements[idx][1]) *
                                        old_solution[0] +
                                        function_derivative_u(x, last_elements[idx][0], last_elements[idx][1]) *
                                        old_solution[1]};
                        }

                );
            } else {
                const auto last_element = ode_solver.get_last_element();
                auto func =
                        [this, last_element]
                                (const vector<solution_t> &old_solution, double x) -> vector<solution_t>
                        {
                            assert(old_solution.size() == 2);

                            return {    old_solution[1],
                                        function_derivative_y( x, last_element[0], last_element[1]) *
                                        old_solution[0] +
                                        function_derivative_u( x, last_element[0], last_element[1]) *
                                        old_solution[1]};
                        };
                ode_solver_for_s_derivative.set_function( func );
                if(have_starting_method){
                    ode_solver_for_s_derivative.starting_method->set_function( func );
                }
            }
            const auto result = ode_solver.get_next_step( step );
            const auto another_result = ode_solver_for_s_derivative.get_next_step( step);
            return {result, another_result};
        }

        [[nodiscard]] bool is_accuracy_reached( double epsilon = 0.000'1) {
            return std::abs(ode_solver.get_last_element()[0] - conditions.right_value) < epsilon;
        }

        solution_t get_new_angle() {
            const auto right_border_result = ode_solver.get_last_element();
            const auto result = current_angle
                                - inverse(ode_solver_for_s_derivative.get_last_element()[0])
                                * ( right_border_result[0] - conditions.right_value );
            return result;
        }

        auto get_solution_error(){
            return ode_solver.get_last_element()[0] - conditions.right_value;
        }

        void restart( const solution_t& angle ) {
            current_angle = angle;

            ode_solver.clear_buffer();
            ode_solver_for_s_derivative.clear_buffer();

            ode_solver.set_initial_conditions({conditions.left_value, current_angle});
            auto func = [this](const vector<solution_t> &old_solution, double x) -> vector<solution_t> {
                assert(old_solution.size() == 2);
                return {old_solution[1], this->function(x, old_solution[0], old_solution[1])};
            };
            ode_solver.set_function(func);
            ode_solver.set_time(conditions.left_border);
            if constexpr (precision_order != 1) {
                ode_solver.starting_method->set_initial_conditions({conditions.left_value, current_angle});
                ode_solver.starting_method->set_function(func);
                ode_solver.starting_method->set_time(conditions.left_border);
            }

            if constexpr (std::is_same_v<solution_t, double>) {
                ode_solver_for_s_derivative.set_initial_conditions({0.0, 1.0});
                if constexpr (precision_order != 1) {
                    ode_solver_for_s_derivative.starting_method->set_initial_conditions({0.0, 1.0});
                }
            } else {
                solution_t first = conditions.left_value;
                for (auto &x : first) {
                    first = 0.0;
                }
                solution_t second = conditions.left_value;
                for (auto &x : second) {
                    second = 1.0;
                }
                ode_solver_for_s_derivative.set_initial_conditions({first, second});
                if constexpr (precision_order != 1) {
                    ode_solver_for_s_derivative.starting_method->set_initial_conditions({first, second});
                }
            }
            ode_solver_for_s_derivative.set_time(conditions.left_border);
            if constexpr (precision_order != 1) {
                ode_solver_for_s_derivative.starting_method->set_time(conditions.left_border);
            }
        }
    };
}
#endif //MAIN_CPP_BOUNDARY_VALUE_PROBLEM_HPP