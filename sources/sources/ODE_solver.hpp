#ifndef MAIN_CPP_ODE_SOLVER_HPP
#define MAIN_CPP_ODE_SOLVER_HPP

#include <functional>
#include <type_traits>

#include "solution_t_concept.hpp"
#include "ring_buffer.hpp"

namespace Numerical_methods{

    template< solution_t_concept solution_t, std::size_t buffer_size>
    class ODE_solver {
    public:
        using function_t = std::function< solution_t(const solution_t&, double)>;
    private:
        virtual return_t<solution_t> get_next_step_impl( double step_, double time_ ) = 0;
    protected:
        ring_buffer<vector<double>, buffer_size> buffer;
        function_t function;
        double time = 0.0;
        solution_t initial_conditions;
    public:
        ODE_solver( const solution_t& conditions, const function_t& function_, double time_ = 0.0 )
                : initial_conditions(conditions),
                  function(function_),
                  time(time_),
                  buffer(ring_buffer<solution_t, buffer_size>())
                {
            buffer.push( initial_conditions );
        };
        [[nodiscard]] double current_time() const noexcept {
            return time;
        }
        return_t<solution_t> get_next_step( double step ) {
            time += step;
            return this->get_next_step_impl( step, time );
        }
    };

    template< solution_t_concept solution_t>
    class one_step_solver : public ODE_solver<solution_t, 1> {
    public:
        using initial_state_t = solution_t;
        using typename ODE_solver<solution_t, 1>::function_t;
        using ODE_solver<solution_t, 1>::ODE_solver;
        using ODE_solver<solution_t, 1>::function;

        return_t<solution_t> make_step (
                const initial_state_t &initial_state,
                double step,
                double time ) const {
            return this->make_step_impl( initial_state, step, time);
        };
    private:
        virtual return_t<solution_t> make_step_impl(
                const initial_state_t &initial_state,
                double step,
                double time ) const = 0;

        return_t<solution_t> get_next_step_impl( double step_, double time_ ){
            return_t<solution_t> result;
            const initial_state_t old_state = ODE_solver<solution_t, 1>::buffer.first();
            result = this->make_step( old_state, step_, time_ );
            const solution_t next_value = result.solution;
            ODE_solver<solution_t, 1>::buffer.push( next_value );
            return result;
        }
    };

    template<solution_t_concept solution_t, std::size_t number_of_steps>
    class multisteps_solver : public ODE_solver<solution_t, number_of_steps - 1>{
    public:
        one_step_solver<solution_t>* starting_method = nullptr;
        using initial_state_t = std::vector<solution_t>;
        using typename ODE_solver<solution_t, number_of_steps - 1>::function_t;
        using ODE_solver<solution_t, number_of_steps - 1>::ODE_solver;
        using ODE_solver<solution_t, number_of_steps - 1>::function;
        void set_starting_method( one_step_solver<solution_t>& method ){
            starting_method = &method;
        };
        return_t<solution_t> make_step (
                const initial_state_t &initial_state,
                double step,
                double time ) const {
            return this->make_step_impl( initial_state, step, time);
        };
        ~multisteps_solver(){
            if( starting_method) {
                //delete starting_method;
                starting_method = nullptr;
            }
        }
    private:
        virtual return_t<solution_t> make_step_impl(
                const initial_state_t &initial_state,
                double step,
                double time ) const = 0;

        return_t<solution_t> get_next_step_impl( double step_, double time_ ){
            return_t<solution_t> result;
            if ( !ODE_solver<solution_t, number_of_steps - 1>::buffer.is_filled() ) {
                const auto old_state = ODE_solver<solution_t, number_of_steps - 1>::buffer.first();
                result = starting_method->make_step( old_state, step_, time_ );
                ODE_solver<solution_t, number_of_steps - 1>::buffer.push( result.solution );
            } else {
                const auto temp = ODE_solver<solution_t, number_of_steps - 1>::buffer.current_state();
                const std::vector<solution_t> old_state{ temp.rbegin(), temp.rend() };
                result = this->make_step( old_state, step_, time_ );
                ODE_solver<solution_t, number_of_steps - 1>::buffer.push( result.solution );
            }
            return result;
        }
    };
    template< solution_t_concept solution_t,
            template <solution_t_concept> class SNAE_method = fixed_point_iterations_method,
                    std::enable_if_t<
                        std::is_base_of_v<  SNAE_solver<solution_t>,
                            fixed_point_iterations_method<solution_t>
                     >, bool> = true>
    struct Implicit{
        static SNAE_method<solution_t> create_SNAE_solver() {
            return SNAE_method<solution_t>();
        }
    };

    struct Explicit{};


}

#endif //MAIN_CPP_ODE_SOLVER_HPP
