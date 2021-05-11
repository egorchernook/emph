#ifndef MAIN_CPP_ONE_STEP_METHODS_HPP
#define MAIN_CPP_ONE_STEP_METHODS_HPP

#include <cstdlib>
#include <limits>
#include <functional>
#include <memory>

#include "solution_t_concept.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "SNAE_solver.hpp"
#include "ODE_solver.hpp"

namespace Numerical_methods{

    template<std::size_t number_of_stages>
    struct Butcher_table {
        vector<double> c_vector;
        vector<double> b_vector;
        matrix<double,number_of_stages> a_matrix;
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
    };
    template<std::size_t number_of_stages>
    struct extended_Butcher_table : private Butcher_table<number_of_stages> {
        vector<double> extra_b_vector;

        using Butcher_table<number_of_stages>::c_vector;
        using Butcher_table<number_of_stages>::b_vector;
        using Butcher_table<number_of_stages>::a_matrix;

        extended_Butcher_table(
                const matrix<double, number_of_stages> &matrix_a,
                const vector<double> &vector_b,
                const vector<double> &vector_c,
                const vector<double> &vector_b2) :
                Butcher_table<number_of_stages>( matrix_a, vector_b, vector_c),
                extra_b_vector(vector_b2)
        {
            assert( vector_c.size() == vector_b.size() );
            assert( matrix_a.height() == matrix_a.width() );
            assert( vector_c.size() == matrix_a.height() );
            assert( vector_b2.size() == vector_b.size() );
        }
    };

    template< solution_t_concept solution_t, std::size_t number_of_stages, std::size_t precision_order>
    class Runge_Kutta_method : public one_step_solver<solution_t> {
    public:
        std::unique_ptr<Butcher_table<number_of_stages>> butcher_table_ptr;
        using typename one_step_solver<solution_t>::function_t;
        using typename one_step_solver<solution_t>::initial_state_t;
        using one_step_solver<solution_t>::one_step_solver;
        using one_step_solver<solution_t>::function;

        void set_Butcher_table( Butcher_table<number_of_stages>& table ) {
            for( int i = 0; i < table.a_matrix.height(); ++i) {
                for( int j = i; j < table.a_matrix.width(); ++j) {
                    if( std::is_base_of_v<Implicit<solution_t>,decltype(*this)>)
                        assert( table.a_matrix[i][j] == 0 );
                }
            }
            butcher_table_ptr = std::make_unique<Butcher_table<number_of_stages>>(table);
        }
    };

    template< solution_t_concept solution_t , std::size_t number_of_stages, std::size_t precision_order>
    class implicit_Runge_Kutta_method : public Runge_Kutta_method< solution_t, number_of_stages, precision_order>,
                                        public Implicit<solution_t> {
    public:
        using base_t = Runge_Kutta_method< solution_t, number_of_stages, precision_order>;
        using base_t::butcher_table_ptr;
        using base_t::Runge_Kutta_method;
        using typename base_t::function_t;
        using initial_state_t = typename base_t::initial_state_t;
        using base_t::function;
    private:
        return_t<solution_t> make_step_impl(
                const initial_state_t &initial_state,
                double step,
                double time ) const override
        {

            vector<solution_t> k_vector( number_of_stages, initial_state);

            auto snae_function = [this, &initial_state, &time, &step]
                    ( const vector<solution_t>& old_solution ) noexcept -> vector<solution_t> {
                vector<solution_t> result = old_solution;
                for (int idx = 0; idx < old_solution.size(); ++idx) {
                    solution_t solution_increments = butcher_table_ptr->a_matrix[idx][0] * old_solution[0];
                    for (int idy = 1; idy < number_of_stages; idy++) {
                        solution_t term = butcher_table_ptr->a_matrix[idx][idy] * old_solution[idy];
                        solution_increments += term;
                    }
                    result[idx] = base_t::function( initial_state + step * solution_increments,
                                                    time + butcher_table_ptr->c_vector[idx] * step );
                }
                return result;
            };

            auto solver = Implicit< vector<solution_t> >::create_SNAE_solver();

            const auto k_vector_result = solver.solve( k_vector, snae_function );

            k_vector = k_vector_result.solution;

            solution_t result = initial_state / step;
            for (int idx = 0; idx < number_of_stages; idx++) {
                result += butcher_table_ptr->b_vector[idx] * k_vector[idx];
            }
            result *= step;

            return { result, k_vector_result.error};
        }
    };

}
#endif //MAIN_CPP_ONE_STEP_METHODS_HPP
