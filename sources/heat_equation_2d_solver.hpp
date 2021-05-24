#ifndef MAIN_CPP_HEAT_EQUATION_2D_SOLVER_HPP
#define MAIN_CPP_HEAT_EQUATION_2D_SOLVER_HPP

#include <functional>
#include <utility>
#include <vector>

#include "solution_t_concept.hpp"
#include "vector.hpp"
#include "matrix.hpp"

namespace Numerical_methods{

    template<solution_t_concept solution_t, std::size_t dimension>
    class boundary_conditions {
    public:
        using coordinate_t = vector<double>;
        using function_t = std::function< solution_t(double, double)>;
    private:
        virtual function_t get_boundary_right_part_function_impl(  const heat_equation_solver& solver ) = 0;
    public:
        const std::vector<function_t> left_values;
        const std::vector<function_t> right_values;
        const coordinate_t left_border;
        const coordinate_t right_border;
        boundary_conditions(
                            std::vector<function_t> left_values_,
                            std::vector<function_t> right_values_,
                            coordinate_t left_border_,
                            coordinate_t right_border_)
                            :
                            left_values(left_values_),
                            right_values(right_values_),
                            left_border(std::move(left_border_)),
                            right_border(std::move(right_border_))
                            {
            assert( left_border.size() == dimension);
            assert( right_border.size() == dimension);
            assert( left_values_.size() == dimension);
            assert( right_values_.size() == dimension);
        }
        function_t get_boundary_right_part_function( const heat_equation_solver& solver) {
            return this->get_boundary_right_part_function_impl( solver );
        };
    };

    template<solution_element_t solution_t>
    class Dirichlet_conditions : public boundary_conditions<solution_t, 1>{
    public:
        using boundary_conditions<solution_t, 1>::boundary_conditions();
    };

    template<solution_element_t solution_t>
    class Neumann_conditions : public boundary_conditions<solution_t, 1>{
    public:
        using boundary_conditions<solution_t, 1>::boundary_conditions();
    };

    template<solution_element_t solution_t>
    class mixed_conditions : public boundary_conditions<solution_t, 1>{
    public:
        using boundary_conditions<solution_t, 1>::boundary_conditions();
    };

    template<solution_element_t solution_t>
    class first_mixed_conditions : public mixed_conditions<solution_t, 1>{
    public:
        using mixed_conditions<solution_t, 1>::mixed_conditions();
    };

    template<solution_element_t solution_t>
    class second_mixed_conditions : public mixed_conditions<solution_t, 1>{
    public:
        using mixed_conditions<solution_t, 1>::mixed_conditions();
    };

    template<solution_element_t solution_t>
    class first_and_second_conditions : public boundary_conditions<solution_t, 2>{
    public:
        using base_t = boundary_conditions<solution_t, 2>;
    private:

    public:
        using boundary_conditions<solution_t, 2>::boundary_conditions();

    };


    class heat_equation_solver{};

    template< template <solution_t_concept = matrix<double>, std::size_t> class solver_t
                = Adams_Moulton_method>
    class heat_equation_2d_solver : public heat_equation_solver {
    public:
        using coordinate_t = vector<double>;
        using solution_t = double;
    private:
        solver_t<> ode_solver;
        std::function<solution_t(const coordinate_t&)> initial_conditions;
        std::function<solution_t(const coordinate_t&)> c_function;
        std::function<solution_t(const coordinate_t&)> lambda_function;
        std::function<solution_t(const coordinate_t&, double)> f_function;
    public:
        const boundary_conditions<double,2> *const conditions = nullptr;
        double time;
        double time_step;
        const coordinate_t lattice_steps_sizes;
        heat_equation_2d_solver(
                                std::function<solution_t(const coordinate_t&)>  initial_conditions_,
                                std::function<solution_t(const coordinate_t&)>  c_function_,
                                std::function<solution_t(const coordinate_t&)>  lambda_function_,
                                std::function<solution_t(const coordinate_t&, double)>  f_function_,
                                const boundary_conditions<double, dimension>& conditions_,
                                const coordinate_t& steps_sizes_,
                                solver_t<>& solver,
                                double initial_time_,
                                double time_step)
                                :   initial_conditions(std::move(initial_conditions_)),
                                    c_function(std::move(c_function_)),
                                    lambda_function(std::move(lambda_function_)),
                                    f_function(std::move(f_function_)),
                                    conditions(conditions_),
                                    lattice_steps_sizes(steps_sizes_),
                                    ode_solver(solver),
                                    time(initial_time_),
                                    time_step(time_step)
                                {
            assert( steps_sizes_.size() == 2 );
            const matrix<solution_t> initial_layer( conditions_.right_border[0] / lattice_steps_sizes[0],
                                                    conditions_.right_border[1] / lattice_steps_sizes[1] );
            for( int i = 0; i < ode_solver.height(); ++i){
                for( int j = 0; j < ode_solver.width(); ++j){
                    initial_layer[i][j] = initial_conditions_( {i * lattice_steps_sizes[0],
                                                                j * lattice_steps_sizes[1]});
                }
            }
            ode_solver.set_initial_conditions(initial_layer);
            ode_solvers[i][j].set_time( time );
                                }

        return_t<matrix<solution_t>> get_next_step() {
            const double h_x =  lattice_steps_sizes[0];
            const double h_y =  lattice_steps_sizes[1];
            ode_solver.set_function(
                    [this, h_x, h_y, time, time_step]
                    ( const matrix<solution_t>& old_solution, double t) -> matrix<solution_t> {
                        auto result = old_solution;
                        for( int i = 1; i < old_solution.height() - 1; ++i) {
                            for (int j = 1; j < old_solution.width() - 1; ++j) {
                                const auto middle_point = old_solution[i][j]
                                result[i][j] =
                                        1.0 / this->c_function( {i * h_x, j * h_y} )
                                            * ( 1.0 / (h_x * h_x)
                                                * ( lambda_function( {i * h_x + 0.5, j * h_y} )
                                                    * ( old_solution[i + 1][j] - middle_point )
                                                - lambda_function( {i * h_x - 0.5, j * h_y} )
                                                    * ( old_solution[i - 1][j] - middle_point ) )
                                            + 1.0 / (h_y * hy_)
                                                * ( lambda_function( {i * h_x, j * h_y + 0.5} )
                                                    * ( old_solution[i][j + 1] - middle_point )
                                                - lambda_function( {i * h_x, j * h_y - 0.5} )
                                                    * ( old_solution[i][j - 1] - middle_point ) ) )
                                        + this->f_function( {i * h_x, j * h_y} );
                                }
                            }

                        }
                    return result;
            );
            time += time_step;
        }

        ~heat_equation_2d_solver() noexcept {
            conditions = nullptr;
        }
    };

}
#endif //MAIN_CPP_HEAT_EQUATION_2D_SOLVER_HPP
