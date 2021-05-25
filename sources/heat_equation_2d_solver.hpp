#ifndef MAIN_CPP_HEAT_EQUATION_2D_SOLVER_HPP
#define MAIN_CPP_HEAT_EQUATION_2D_SOLVER_HPP

#include <functional>
#include <utility>
#include <vector>

#include "solution_t_concept.hpp"
#include "vector.hpp"
#include "matrix.hpp"

namespace Numerical_methods{


    class heat_equation_solver{
    public:
        using coordinate_t = vector<double>;
        using solution_t = double;
        using function_t = std::function<solution_t(const coordinate_t&)>;
        using solution_function_t = std::function<solution_t(const coordinate_t&, double)>;
    };

    template< std::size_t precision_order,
                template <solution_t_concept = matrix<double>, std::size_t = precision_order> class solver_t
                = Adams_Moulton_method>
    class heat_equation_2d_solver;

    template<solution_element_t solution_t, std::size_t dimension>
    class boundary_conditions {
    public:
        using coordinate_t = vector<double>;
        using function_t = std::function< solution_t(double, double)>;
    private:
        [[nodiscard]] virtual std::function< matrix<heat_equation_solver::solution_t>
                                        ( const matrix<heat_equation_solver::solution_t>&, double)>
                                        get_function_impl (
                                                const heat_equation_solver::function_t& c_function,
                                                const heat_equation_solver::function_t& lambda_function,
                                                const heat_equation_solver::solution_function_t& f_function,
                                                std::size_t height,
                                                std::size_t width ) const= 0;
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
        [[nodiscard]] std::function< matrix<heat_equation_solver::solution_t>
                ( const matrix<heat_equation_solver::solution_t>&, double)>
            get_function( const heat_equation_solver::function_t& c_function,
                                           const heat_equation_solver::function_t& lambda_function,
                                           const heat_equation_solver::solution_function_t& f_function,
                                           std::size_t height,
                                           std::size_t width ) const {
                return this->get_function_impl( c_function, lambda_function, f_function, height, width );
        };
    };

    template<solution_element_t solution_t>
    class Dirichlet_conditions : public boundary_conditions<solution_t, 1>{
    public:
        using boundary_conditions<solution_t, 1>::boundary_conditions;
    };

    template<solution_element_t solution_t>
    class Neumann_conditions : public boundary_conditions<solution_t, 1>{
    public:
        using boundary_conditions<solution_t, 1>::boundary_conditions;
    };

    template<solution_element_t solution_t>
    class mixed_conditions : public boundary_conditions<solution_t, 1>{
    public:
        using boundary_conditions<solution_t, 1>::boundary_conditions;
    };

    template<solution_element_t solution_t>
    class first_mixed_conditions : public mixed_conditions<solution_t>{
    public:
        using mixed_conditions<solution_t>::mixed_conditions;
    };

    template<solution_element_t solution_t>
    class second_mixed_conditions : public mixed_conditions<solution_t>{
    public:
        using mixed_conditions<solution_t>::mixed_conditions;
    };

    template<solution_element_t solution_t>
    class first_and_second_conditions : public boundary_conditions<solution_t, 2>{
    public:
        using base_t = boundary_conditions<solution_t, 2>;
    private:
        [[nodiscard]] std::function< matrix<heat_equation_solver::solution_t>
                ( const matrix<heat_equation_solver::solution_t>&, double)>
                get_function_impl (
                                        const heat_equation_solver::function_t& c_function,
                                        const heat_equation_solver::function_t& lambda_function,
                                        const heat_equation_solver::solution_function_t& f_function,
                                        std::size_t height,
                                        std::size_t width ) const override {

            const auto h_x = base_t::right_border[0] - base_t::left_border[0];
            const auto h_y = base_t::right_border[1] - base_t::left_border[1];

            const auto left_x = [this]( double y, double time) -> double {
                return base_t::left_values[0](y, time);
            };
            const auto right_x = [this]( double y, double time) -> double {
                return base_t::right_values[0](y, time);
            };
            const auto Mu_b = base_t::left_values[1];
            const auto Mu_u = base_t::right_values[1];
            const auto result = [&c_function,&lambda_function,&f_function,h_x,h_y,left_x,right_x,Mu_b,Mu_u]
                    (const matrix<heat_equation_solver::solution_t> &old_solution, double time)
                    -> matrix<heat_equation_solver::solution_t> {
                auto result = old_solution;
                for (int i = 1; i < old_solution.height() - 1; ++i) {
                    for (int j = 1; j < old_solution.width() - 1; ++j) {
                        const auto middle_point = old_solution[i][j];
                        result[i][j] =
                                1.0 / c_function({i * h_x, j * h_y})
                                * (1.0 / (h_x * h_x)
                                   * (lambda_function({i * h_x + 0.5, j * h_y})
                                      * (old_solution[i + 1][j] - middle_point)
                                      - lambda_function({i * h_x - 0.5, j * h_y})
                                        * (middle_point - old_solution[i - 1][j] ))
                                   + 1.0 / (h_y * h_y)
                                     * (lambda_function({i * h_x, j * h_y + 0.5})
                                        * (old_solution[i][j + 1] - middle_point)
                                        - lambda_function({i * h_x, j * h_y - 0.5})
                                          * (middle_point - old_solution[i][j - 1] )))
                                + f_function({i * h_x, j * h_y}, time);
                    }
                }
                for( int i = 0; i < result.width(); ++i ) {
                    result[0][i] = left_x( i * h_y , time );
                    result[result.height() - 1][i] = right_x( i * h_y , time );
                }
                for( int i = 1; i < result.height() - 1; ++i ) {
                    int j = 0;
                    auto middle_point = old_solution[i][j];
                    result[i][j] =
                            1.0 / c_function({i * h_x, j * h_y})
                            * (1.0 / (h_x * h_x)
                               * (lambda_function({i * h_x + 0.5, j * h_y})
                                  * (old_solution[i + 1][j] - middle_point)
                                  - lambda_function({i * h_x - 0.5, j * h_y})
                                    * (middle_point - old_solution[i - 1][j] ))
                               + 1.0 / (h_y * h_y)
                                 * (lambda_function({i * h_x, j * h_y + 0.5})
                                    * (old_solution[i][j + 1] - middle_point)
                                    - lambda_function({i * h_x, j * h_y - 0.5})
                                      * (middle_point
                                        - ( old_solution[i][j + 1] - 2 * h_y * Mu_b( i * h_x, time ) ) ) ) )
                            + f_function({i * h_x, j * h_y}, time);

                    j = static_cast<int>(result.width()) - 1;
                    middle_point = old_solution[i][j];
                    result[i][j] =
                            1.0 / c_function({i * h_x, j * h_y})
                            * (1.0 / (h_x * h_x)
                               * (lambda_function({i * h_x + 0.5, j * h_y})
                                  * (old_solution[i + 1][j] - middle_point)
                                  - lambda_function({i * h_x - 0.5, j * h_y})
                                    * (middle_point - old_solution[i - 1][j] ))
                               + 1.0 / (h_y * h_y)
                                 * (lambda_function({i * h_x, j * h_y + 0.5})
                                    * ( ( old_solution[i][j - 1] + 2 * h_y * Mu_u( i * h_x,time) )
                                        - middle_point)
                                    - lambda_function({i * h_x, j * h_y - 0.5})
                                      * (middle_point - old_solution[i][j - 1] )))
                            + f_function({i * h_x, j * h_y}, time);
                }
                return result;
            };
            return result;
        };
    public:
        using boundary_conditions<solution_t, 2>::boundary_conditions;
    };

    template< std::size_t precision_order,
                template <solution_t_concept = matrix<double>,
    std::size_t = precision_order> class solver_t>
    class heat_equation_2d_solver : public heat_equation_solver {
    public:
        using base_t = heat_equation_solver;
        using base_t::coordinate_t;
        using base_t::solution_t;
        using base_t::function_t;
        using base_t::solution_function_t;
    private:
        solver_t<> ode_solver;
        function_t initial_conditions;
        function_t c_function;
        function_t lambda_function;
        solution_function_t f_function;
    public:
        const boundary_conditions<double, 2> *const conditions = nullptr;
        double time;
        double time_step;
        const coordinate_t lattice_steps_sizes;
        heat_equation_2d_solver(
                                function_t  initial_conditions_,
                                function_t  c_function_,
                                function_t  lambda_function_,
                                solution_function_t  f_function_,
                                const boundary_conditions<double, 2>* const& conditions_,
                                const coordinate_t& steps_sizes_,
                                solver_t<>& solver,
                                double initial_time_,
                                double time_step)
                                :   initial_conditions( initial_conditions_),
                                    c_function( c_function_),
                                    lambda_function( lambda_function_),
                                    f_function( f_function_),
                                    conditions(conditions_),
                                    lattice_steps_sizes(steps_sizes_),
                                    ode_solver(solver),
                                    time(initial_time_),
                                    time_step(time_step)
                                {
            assert( steps_sizes_.size() == 2 );

            const matrix<solution_t> initial_layer( conditions_->right_border[0] / lattice_steps_sizes[0],
                                                    conditions_->right_border[1] / lattice_steps_sizes[1] );
            for( int i = 0; i < initial_layer.height(); ++i){
                for( int j = 0; j < initial_layer.width(); ++j){
                    initial_layer[i][j] = initial_conditions_( {i * steps_sizes_[0],
                                                                j * steps_sizes_[1]});
                }
            }
            ode_solver.set_initial_conditions(initial_layer);
            ode_solver.set_time( time );
            ode_solver.set_function( conditions->get_function(
                    this->c_function,
                    this->lambda_function,
                    this->f_function,
                    initial_layer.height(),
                    initial_layer.width())
                    );
                                }

        return_t<matrix<solution_t>> get_next_step() {
            const auto result = ode_solver.get_next_step( time_step );
            time += time_step;
            return result;
        }
    };

}
#endif //MAIN_CPP_HEAT_EQUATION_2D_SOLVER_HPP
