#ifndef MAIN_CPP_SNAE_SOLVER_HPP
#define MAIN_CPP_SNAE_SOLVER_HPP

#include <cstdlib>
#include <functional>
#include <type_traits>

#include "solution_t_concept.hpp"
#include "min_and_max_functions.hpp"
#include "matrix.hpp"
#include "expand_function.hpp"

namespace Numerical_methods {


    template<solution_element_t solution_t>
    class SNAE_solver {
    public:
        using function_t = std::function<solution_t(const solution_t &)>;
    private:
        virtual return_t<solution_t> solve_impl(const solution_t &) = 0;
    protected:
        function_t function;
    public:
        double methods_error = 0.000'1;

        SNAE_solver(  const function_t& function_, const double &error = 0.000'1)
                    : methods_error(error),
                      function(function_){}

        return_t<solution_t> solve(const solution_t &initial_approximation) {
            return this->solve_impl(initial_approximation);
        }
    };

    template<solution_t_concept solution_t>
    class fixed_point_iterations_method : public SNAE_solver<solution_t> {
    public:
        using base_t = SNAE_solver<solution_t>;
        using typename base_t::function_t;
        using SNAE_solver<solution_t>::SNAE_solver;
    private:
        return_t<solution_t> solve_impl( const solution_t &initial_state) override {
            solution_t old_solution = initial_state;
            solution_t new_solution = initial_state;

            auto get_local_error = [&new_solution, &old_solution]() noexcept -> double {
                return Numerical_methods::max(Numerical_methods::abs(new_solution - old_solution));
            };

            do {
                old_solution = new_solution;
                new_solution = base_t::function( old_solution );
            } while (get_local_error() - this->methods_error > 0.0);

            return {new_solution, get_local_error()};
        }
    };

    template<solution_t_concept solution_t>
    class Seidels_method : public SNAE_solver<solution_t> {
    public:
        using base_t = SNAE_solver<solution_t>;
        using typename base_t::function_t;
        using SNAE_solver<solution_t>::SNAE_solver;
    private:
        return_t<solution_t> solve_impl(
                const solution_t &initial_state) override {
            solution_t old_solution = initial_state;
            solution_t new_solution = initial_state;

            auto get_local_error = [&new_solution, &old_solution]() noexcept -> double {
                return Numerical_methods::max(Numerical_methods::abs(new_solution - old_solution));
            };

            do {
                old_solution = new_solution;
                for (int idx = 0; idx < old_solution.size(); ++idx) {
                    new_solution[idx] = base_t::function(new_solution)[idx];
                }
            } while (get_local_error() > this->methods_error);

            return {new_solution, get_local_error()};
        }
    };

    template<solution_t_concept solution_t>
    class Newtons_method : public SNAE_solver<solution_t>{
    public:
        using base_t = SNAE_solver<solution_t>;
        using typename SNAE_solver<solution_t>::function_t;
        using value_t = std::decay_t<decltype(std::declval<solution_t>()[0])>;
    private:
        std::vector<function_t> function_derivatives;
    public:
        Newtons_method( const function_t& function_, const std::vector<function_t>& function_derivatives_, const double &error = 0.000'1)
            : base_t(function_, error),
            function_derivatives(function_derivatives_){}
    private:
        return_t<solution_t> solve_impl(
                const solution_t & initial_approximation ) override {
            solution_t error{};

            const int solution_size = size(initial_approximation);

            solution_t old_solution = initial_approximation;
            solution_t new_solution = initial_approximation;

            auto get_local_error = [&new_solution, &old_solution]() noexcept -> double {
                return Numerical_methods::max(Numerical_methods::abs(new_solution - old_solution));
            };

            matrix<value_t> derivative( solution_size);

            do {
                old_solution = new_solution;

                for( int i = 0; i < solution_size; ++i){
                    const auto row = function_derivatives[i](old_solution);
                    for( int j = 0; j < solution_size; ++j) {
                        derivative[i][j] = row[j];
                    }
                }
                const solution_t free_part = -1.0 * base_t::function( old_solution );
                error = Gauss_method(
                            derivative,
                            free_part);

                new_solution = old_solution + error;
            } while ( get_local_error() > this->methods_error );

            return { new_solution, max(error) };
        }
    };
}
#endif //MAIN_CPP_SNAE_SOLVER_HPP
