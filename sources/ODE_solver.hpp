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
        void set_initial_conditions_impl(const solution_t &initial_conditions_) {
            initial_conditions = initial_conditions_;
            buffer.push( initial_conditions );
        }
        void set_function_impl( const function_t& function_ ) {
            function = function_;
        }

        void set_time_impl( double time_ ) {
            time = time_;
        }
        ring_buffer<solution_t, buffer_size> buffer;
        function_t function;
        double time = 0.0;
        solution_t initial_conditions;
    public:
        auto get_buffer_state()  {
            return buffer.current_state();
        }
        void clear_buffer(){
            buffer.clear();
        }
        [[nodiscard]] bool buffer_is_filled()  {
            return buffer.is_filled();
        }
        auto get_last_element()  {
            return buffer.first();
        }
        ODE_solver() : buffer(ring_buffer<solution_t, buffer_size>()) {}

        virtual void set_initial_conditions(const solution_t &initial_conditions_) {
            set_initial_conditions_impl( initial_conditions_);
        }
        virtual void set_function( const function_t& function_ ) {
            set_function_impl( function_);
        }

        virtual void set_time( double time_ ) {
            set_time_impl( time_);
        }
        [[nodiscard]] double get_time() const noexcept {
            return time;
        }

        ODE_solver( const solution_t& conditions, const function_t& function_, double time_ = 0.0 )
                : initial_conditions(conditions),
                  function(function_),
                  time(time_),
                  buffer(ring_buffer<solution_t, buffer_size>())
                {
            buffer.push( initial_conditions );
        }

        return_t<solution_t> get_next_step( double step ) {
            time += step;
            return this->get_next_step_impl( step, time );
        }
        [[nodiscard]] int size() const {
            return buffer_size;
        }
    };

    template< solution_t_concept solution_t>
    class one_step_solver : public ODE_solver<solution_t, 1> {
    public:
        using base_t = ODE_solver<solution_t, 1>;
        using initial_state_t = solution_t;
        using typename base_t::function_t;
        using ODE_solver<solution_t, 1>::ODE_solver;
        using base_t::function;

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
            const initial_state_t old_state = base_t::buffer.first();
            result = this->make_step( old_state, step_, time_ );
            const solution_t next_value = result.solution;
            base_t::buffer.push( next_value );
            return result;
        }

    };

    template<solution_t_concept solution_t, std::size_t number_of_steps>
    class multisteps_solver : public ODE_solver<solution_t, number_of_steps - 1>{
    public:
        using base_t = ODE_solver<solution_t, number_of_steps - 1>;
        one_step_solver<solution_t>* starting_method = nullptr;
        using initial_state_t = std::vector<solution_t>;
        using typename base_t::function_t;
        using ODE_solver<solution_t, number_of_steps - 1>::ODE_solver;
        using base_t::function;
        void set_starting_method( one_step_solver<solution_t>& method ){
            starting_method = &method;
        };
        return_t<solution_t> make_step (
                const initial_state_t &initial_state,
                double step,
                double time ) const {
            return this->make_step_impl( initial_state, step, time);
        };
        void set_initial_conditions(const solution_t &initial_conditions_) override {
            if( starting_method){
                starting_method->set_initial_conditions( initial_conditions_);
            }
            base_t::set_initial_conditions_impl( initial_conditions_);
        }
        void set_function( const function_t& function_ ) override {
            if( starting_method){
                starting_method->set_function( function_);
            }
            base_t::set_function_impl( function_);
        }
        void set_time( double time_ ) override {
            if( starting_method){
                starting_method->set_time( time_);
            }
            base_t::set_time_impl( time_);
        }
        ~multisteps_solver(){
            if( starting_method) {
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
            if ( !base_t::buffer.is_filled() ) {
                const auto old_state = base_t::buffer.first();
                result = starting_method->make_step( old_state, step_, time_ );
                base_t::buffer.push( result.solution );
            } else {
                const auto temp = base_t::buffer.current_state();
                const std::vector<solution_t> old_state{ temp.rbegin(), temp.rend() };
                result = this->make_step( old_state, step_, time_ );
                base_t::buffer.push( result.solution );
            }
            return result;
        }
    };
    template< solution_t_concept solution_t,
            template <solution_t_concept> class SNAE_method = Seidels_method,
                    std::enable_if_t<
                        std::is_base_of_v<  SNAE_solver<solution_t>,
                            SNAE_method<solution_t>
                     >, bool> = true>
    struct Implicit{
        static SNAE_method<solution_t> create_SNAE_solver( const typename SNAE_method<solution_t>::function_t& function) {
            return SNAE_method<solution_t>( function );
        }
    };

    struct Explicit{};


}

#endif //MAIN_CPP_ODE_SOLVER_HPP
