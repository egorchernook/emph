#ifndef MAIN_CPP_MULTISTEPS_METHODS_HPP
#define MAIN_CPP_MULTISTEPS_METHODS_HPP

#include <cstdlib>
#include <vector>
#include <type_traits>
#include "SNAE_solver.hpp"
#include "solution_t_concept.hpp"
#include "ODE_solver.hpp"

namespace Numerical_methods{

    template<solution_t_concept solution_t, std::size_t number_of_steps>
    class Adams_method : public multisteps_solver<solution_t,number_of_steps> {
    public:
        using typename multisteps_solver<solution_t,number_of_steps>::function_t;
        using typename multisteps_solver<solution_t,number_of_steps>::initial_state_t;
        using multisteps_solver<solution_t,number_of_steps>::multisteps_solver;
        using multisteps_solver<solution_t,number_of_steps>::function;
        std::vector<double> coefficients{0.0};
        void set_coefficients( const std::vector<double> &coefficients_ ) {
            assert( coefficients_.size() == number_of_steps);
            coefficients = coefficients_;
        }
    };

    template< solution_t_concept solution_t, std::size_t number_of_steps>
    class Adams_Moulton_method : public Adams_method<solution_t,number_of_steps>,
                                 public Implicit<solution_t> {
    public:
        using base_t = Adams_method<solution_t,number_of_steps>;
        using base_t::Adams_method;
        using typename base_t::function_t;
        using typename base_t::initial_state_t;
        using base_t::function;
    private:
        return_t<solution_t> make_step_impl(
                const initial_state_t &initial_state,
                double step,
                double time ) const override {

            assert(initial_state.size() == number_of_steps - 1);

            solution_t constant_part = this->coefficients[1] * function(
                                        initial_state[0],
                                        time - step );
            for (int idx = 2; idx < number_of_steps; ++idx ) {
                const solution_t temp = this->coefficients[idx] * function(
                        initial_state[ idx - 1 ],
                        time - step * idx);
                constant_part += temp;
            }

            constant_part *= step;
            constant_part += initial_state[0];

            const double beta_zero = this->coefficients[0];
            auto right_part_function = [&]( const solution_t& solution) noexcept -> solution_t {
                        return step * beta_zero * function( solution, time) + constant_part;
            };

            auto solver = Implicit<solution_t>::create_SNAE_solver();

            return_t<solution_t> result = solver.solve( initial_state[0] , right_part_function );
            return result;
        }
    };
}
#endif //MAIN_CPP_MULTISTEPS_METHODS_HPP
