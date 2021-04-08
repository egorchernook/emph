#ifndef MAIN_CPP_MULTISTEPS_METHODS_HPP
#define MAIN_CPP_MULTISTEPS_METHODS_HPP

#include <cstdlib>
#include <array>
#include <vector>
#include "SNAE_solver.hpp"
#include "solution_t_concept.hpp"
#include "ODE_solver.hpp"

namespace Numerical_methods{

    template<solution_t_concept solution_t>
    class multisteps_methods : public Cauchy_problem_solver<solution_t> {
    public:
        using initial_state_t = std::vector<solution_t>;
    };

    template<solution_t_concept solution_t, std::size_t precision_order>
    class Adams_method : public multisteps_methods<solution_t> {
    public:
        using function_t = std::function< solution_t(const solution_t&, double)>;
        using typename multisteps_methods<solution_t>::initial_state_t;
    private:
        virtual return_t<solution_t> make_step_impl(const initial_state_t &initial_state,
                                                    const function_t& initial_function,
                                                    double step,
                                                    double time ) const = 0;
    public:
        std::array<double, precision_order> coefficients{0.0};
        Adams_method(const std::array<double, precision_order> &_coefficients) :
            coefficients(_coefficients) {}
        return_t<solution_t> make_step(
                const initial_state_t &initial_state,
                const function_t& initial_function,
                double step,
                double time ) const {
            return this->make_step_impl( initial_state, initial_function, step, time);
        };
    };

    template<solution_t_concept solution_t, std::size_t precision_order>
    class Adams_Moulton_method : public Adams_method<solution_t,precision_order> {
        using Adams_method<solution_t,precision_order>::Adams_method;
        using typename Adams_method<solution_t,precision_order>::function_t;
        using typename Adams_method<solution_t,precision_order>::initial_state_t;
    private:
        return_t<solution_t> make_step_impl(
                const initial_state_t &initial_state,
                const function_t& initial_function,
                double step,
                double time ) const override {

            assert(initial_state.size() == precision_order - 1);

            solution_t constant_part = this->coefficients[1] * initial_function(
                                        initial_state[0],
                                        time - step );
            for (int idx = 2; idx < precision_order; ++idx ) {
                const solution_t temp = this->coefficients[idx] * initial_function(
                                                                                    initial_state[ idx - 1 ],
                                                                                    time - step * idx);
                constant_part += temp;
            }

            constant_part *= step;
            constant_part += initial_state[0];

            const double beta_zero = this->coefficients[0];
            auto right_part_function = [&]( const solution_t& solution) noexcept -> solution_t {
                        return step * beta_zero * initial_function( solution, time) + constant_part;
            };

            fixed_point_iterations_method< solution_t> solver{};

            return_t<solution_t> result = solver.solve( initial_state[0] , right_part_function );
            return result;
        }
    };
}
#endif //MAIN_CPP_MULTISTEPS_METHODS_HPP
