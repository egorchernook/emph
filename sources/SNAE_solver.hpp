#ifndef MAIN_CPP_SNAE_SOLVER_HPP
#define MAIN_CPP_SNAE_SOLVER_HPP

#include <cstdlib>
#include <functional>
#include "solution_t_concept.hpp"
#include "min_and_max_functions.hpp"

namespace Numerical_methods {


    template<solution_t_concept solution_t>
    class SNAE_solver {
    public:
        using function_t = std::function<solution_t(const solution_t &)>;
    private:
        virtual return_t<solution_t> solve_impl(const solution_t &, const function_t &) = 0;

    public:
        double methods_error = 0.000'001;

        SNAE_solver(const double &error = 0.000'001) : methods_error(error) {}

        return_t<solution_t> solve(const solution_t &solution, const function_t &function) {
            return this->solve_impl(solution, function);
        }
    };

    /*
    template<solution_t_concept solution_t>
    class Newtons_method : public SNAE_solving_methods{
    private:
        [[nodiscard]] bool is_local_error_small_enough( double local_error ) const {
            return method_error * 100'000 - local_error * 100'000 > 0.0;
        }
    public:
        using SNAE_solving_methods::SNAE_solving_methods;

        return_t<solution_t> make_step(
                const solution_t & derivative,
                const solution_t & free_part,
                const solution_t & initial_approximation
                ) {

            solution_t error{};

            error = Gauss_method<solution_t, size(free_part)>(
                                                                    derivative,
                                                                    free_part);

            const solution_t new_solution = initial_approximation + error;

            return { new_solution, max(error) };
        }
    };
*/
    template<solution_t_concept solution_t>
    class fixed_point_iterations_method : public SNAE_solver<solution_t> {
    public:
        using typename SNAE_solver<solution_t>::function_t;
        using SNAE_solver<solution_t>::SNAE_solver;
    private:
        return_t<solution_t> solve_impl(
                const solution_t &initial_state,
                const function_t &right_part_function) override {
            solution_t old_solution = initial_state;
            solution_t new_solution = initial_state;

            auto get_local_error = [&new_solution, &old_solution]() noexcept -> double {
                return Numerical_methods::max(Numerical_methods::abs(new_solution - old_solution));
            };

            do {
                old_solution = new_solution;
                new_solution = right_part_function(old_solution);
            } while (get_local_error() - this->methods_error > 0.0);

            return {new_solution, get_local_error()};
        }
    };

    template<solution_t_concept solution_t>
    class Seidels_method : public SNAE_solver<solution_t> {
    public:
        using typename SNAE_solver<solution_t>::function_t;
        using SNAE_solver<solution_t>::SNAE_solver;
    private:
        return_t<solution_t> solve_impl(
                const solution_t &initial_state,
                const function_t &right_part_function) override {
            solution_t old_solution = initial_state;
            solution_t new_solution = initial_state;

            auto get_local_error = [&new_solution, &old_solution]() noexcept -> double {
                return Numerical_methods::max(Numerical_methods::abs(new_solution - old_solution));
            };

            do {
                old_solution = new_solution;
                for (int idx = 0; idx < old_solution.size(); ++idx) {
                    new_solution[idx] = right_part_function(new_solution)[idx];
                }
            } while (get_local_error() > this->methods_error);

            return {new_solution, get_local_error()};
        }
    };
}
#endif //MAIN_CPP_SNAE_SOLVER_HPP
