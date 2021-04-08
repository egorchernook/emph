#ifndef MAIN_CPP_ONE_STEP_METHODS_HPP
#define MAIN_CPP_ONE_STEP_METHODS_HPP

#include <cstdlib>
#include <limits>
#include <functional>

#include "solution_t_concept.hpp"
#include "vector.hpp"
#include "SNAE_solver.hpp"
#include "ODE_solver.hpp"

namespace Numerical_methods{

    template< solution_t_concept solution_t>
    class one_step_ODE_solver : public Cauchy_problem_solver<solution_t> {
    public:
        using initial_state_t = solution_t;
    };

    template<std::size_t number_of_stages>
    struct Butcher_table {
        //vector<double>* extra_b_vector_ptr = nullptr;
        vector<double> c_vector;
        vector<double> b_vector;
        matrix<double, number_of_stages> a_matrix;

        Butcher_table(
                const matrix<double, number_of_stages> &matrix_a,
                const vector<double> &vector_b,
                const vector<double> &vector_c ) :
                    c_vector(vector_c),
                    b_vector(vector_b),
                    a_matrix(matrix_a)
                {
            assert( vector_c.size() == vector_b.size() );
            assert( matrix_a.height() == matrix_a.width() );
            assert( vector_c.size() == matrix_a.height() );
        }
        /*
        Butcher_table(
                const matrix<double, number_of_stages> &matrix_a,
                const vector<double> &vector_b,
                const vector<double> &vector_c,
                vector<double> &ptr_vector_b_extra )
                : extra_b_vector_ptr( &ptr_vector_b_extra )
        {

            assert( vector_c.size() == vector_b.size() );
            assert( matrix_a.height() == matrix_a.width() );
            assert( vector_c.size() == matrix_a.height() );
            assert( extra_b_vector_ptr->size() == vector_b.size() );

            c_vector = vector_c;
            b_vector = vector_b;
            a_matrix = matrix_a;
        }
         */
    };

    template< solution_t_concept solution_t, std::size_t number_of_stages>
    class Runge_Kutta_method : public one_step_ODE_solver<solution_t> {
    public:
        using function_t = std::function< solution_t(solution_t, double)>;
        using initial_state_t = typename one_step_ODE_solver<solution_t>::initial_state_t;
    private:
        virtual return_t<solution_t> make_step_impl(const initial_state_t &initial_state,
                                                    const function_t& initial_function,
                                                    double step,
                                                    double time ) const = 0;
    public:
        return_t<solution_t> make_step(
                const initial_state_t &initial_state,
                const function_t &function,
                double step,
                double time ) const {
            return this->make_step_impl( initial_state, function, step, time);
        }
    };

    template< solution_t_concept solution_t , std::size_t number_of_stages>
    class implicit_Runge_Kutta_method
            : public Runge_Kutta_method< solution_t, number_of_stages>,
              private Butcher_table<number_of_stages> {
    public:
        using typename Runge_Kutta_method< solution_t, number_of_stages>::function_t;
        using initial_state_t = typename Runge_Kutta_method<solution_t,number_of_stages>::initial_state_t;
    private:
        return_t<solution_t> make_step_impl(
                const initial_state_t &initial_state,
                const function_t &function,
                double step,
                double time ) const override
        {

            vector<solution_t> k_vector( number_of_stages, initial_state);
            fixed_point_iterations_method< vector<solution_t>> solver{};

            auto snae_function = [&, this]( const vector<solution_t>& old_solution ) noexcept -> vector<solution_t> {
                vector<solution_t> result = old_solution;
                for (int idx = 0; idx < old_solution.size(); ++idx) {
                    solution_t solution_increments = this->a_matrix[idx][0] * old_solution[0];
                    for (int idy = 1; idy < number_of_stages; idy++) {
                        solution_t term = this->a_matrix[idx][idy] * old_solution[idy];
                        solution_increments += term;
                    }
                    result[idx] = function( initial_state + step * solution_increments,
                                            time + this->c_vector[idx] * step );
                }
                return result;
            };
            return_t< vector<solution_t> > k_vector_result = solver.solve( k_vector, snae_function );

            k_vector = k_vector_result.solution;

            solution_t result = initial_state / step;
            for (int idx = 0; idx < number_of_stages; idx++) {
                result += this->b_vector[idx] * k_vector[idx];
            }
            result *= step;

            return { result, k_vector_result.error};
        }
    public:
        using Butcher_table<number_of_stages>::Butcher_table;
    };

}
#endif //MAIN_CPP_ONE_STEP_METHODS_HPP
