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
        using function_t = std::function<std::decay<solution_t>::type(const coordinate_t&)>;
        using solution_function_t = std::function<std::decay<solution_t>::type(const coordinate_t&, double)>;
    };

    template< std::size_t precision_order,
                template <solution_t_concept = matrix<double>, std::size_t = precision_order> class solver_t
                = Adams_Moulton_method>
    class heat_equation_2d_solver;

    template<solution_element_t solution_t, std::size_t dimension>
    class boundary_conditions {
    private:
        [[nodiscard]] virtual bool is_left_derivative_impl() const = 0;
        [[nodiscard]] virtual bool is_right_derivative_impl() const = 0;
    public:
        using function_t = std::function<
                typename std::decay<solution_t>::type
                        ( const vector<double>&, double)>;
    public:
        const function_t left_function;
        const function_t right_function;
        const double left_border;
        const double right_border;

        [[nodiscard]] bool is_left_derivative() const {
            return this->is_left_derivative_impl();
        };
        [[nodiscard]] bool is_right_derivative() const {
            return this->is_right_derivative_impl();
        };

        boundary_conditions(
                            function_t left_values_,
                            function_t right_values_,
                            double left_border_,
                            double right_border_)
                            :
                            left_function(left_values_),
                            right_function(right_values_),
                            left_border( left_border_),
                            right_border(right_border_)
                            {
        }
    };

    template<solution_element_t solution_t, std::size_t dimension>
    class Dirichlet_conditions : public boundary_conditions<solution_t, dimension>{
    private:
        [[nodiscard]] bool is_left_derivative_impl() const {
            return false;
        }
        [[nodiscard]] bool is_right_derivative_impl() const {
            return false;
        }
    public:
        using boundary_conditions<solution_t, dimension>::boundary_conditions;
    };

    template<solution_element_t solution_t, std::size_t dimension>
    class Neumann_conditions : public boundary_conditions<solution_t, dimension>{
    private:
        [[nodiscard]] bool is_left_derivative_impl() const {
            return true;
        }
        [[nodiscard]] bool is_right_derivative_impl() const {
            return true;
        }
    public:
        using boundary_conditions<solution_t, dimension>::boundary_conditions;
    };

    template<solution_element_t solution_t, std::size_t dimension>
    class mixed_conditions : public boundary_conditions<solution_t, dimension>{
    public:
        using boundary_conditions<solution_t, dimension>::boundary_conditions;
    };

    template<solution_element_t solution_t, std::size_t dimension>
    class first_mixed_conditions : public mixed_conditions<solution_t, dimension>{
    private:
        [[nodiscard]] bool is_left_derivative_impl() const {
            return false;
        }
        [[nodiscard]] bool is_right_derivative_impl() const {
            return true;
        }
    public:
        using mixed_conditions<solution_t, dimension>::mixed_conditions;
    };

    template<solution_element_t solution_t, std::size_t dimension>
    class second_mixed_conditions : public mixed_conditions<solution_t, dimension>{
    private:
        [[nodiscard]] bool is_left_derivative_impl() const {
            return true;
        }
        [[nodiscard]] bool is_right_derivative_impl() const {
            return false;
        }
    public:
        using mixed_conditions<solution_t, dimension>::mixed_conditions;
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
        const boundary_conditions<double, 2> *const x_conditions = nullptr;
        const boundary_conditions<double, 2> *const y_conditions = nullptr;
        double time;
        double time_step;
        const coordinate_t grid_steps_sizes;

        heat_equation_2d_solver(
                function_t initial_conditions_,
                function_t c_function_,
                function_t lambda_function_,
                solution_function_t f_function_,
                const boundary_conditions<double, 2> *const &x_conditions_,
                const boundary_conditions<double, 2> *const &y_conditions_,
                const coordinate_t &steps_sizes_,
                solver_t<> &solver,
                double initial_time_,
                double time_step)
                : initial_conditions(initial_conditions_),
                  c_function(c_function_),
                  lambda_function(lambda_function_),
                  f_function(f_function_),
                  x_conditions(x_conditions_),
                  y_conditions(y_conditions_),
                  grid_steps_sizes(steps_sizes_),
                  ode_solver(solver),
                  time(initial_time_),
                  time_step(time_step) {
            assert(steps_sizes_.size() == 2);

            const std::size_t width =
                    static_cast<int>((y_conditions_->right_border - y_conditions_->left_border) / steps_sizes_[1]) + 1;
            const std::size_t height =
                    static_cast<int>((x_conditions_->right_border - x_conditions_->left_border) / steps_sizes_[0]) + 1;
            const matrix<solution_t> initial_layer(height, width);

            for (int i = 0; i < initial_layer.height(); ++i) {
                for (int j = 0; j < initial_layer.width(); ++j) {
                    const double x = x_conditions_->left_border + i * steps_sizes_[0];
                    const double y = y_conditions_->left_border + j * steps_sizes_[1];
                    initial_layer[i][j] = initial_conditions_({x, y});
                }
            }
            ode_solver.set_initial_conditions(initial_layer);
            ode_solver.set_time(time);

            const auto h_x = steps_sizes_[0];
            const auto h_y = steps_sizes_[1];


            std::vector< std::function< std::decay<double>::type
                    (const solution_t&,
                     const solution_t&,
                     const solution_t&,
                     const solution_t&,
                     const solution_t&,
                     double,
                     double)>> functions( 4);
            if (!x_conditions->is_left_derivative()) {

                functions[0] = [this]( const solution_t& left, const solution_t& right,
                                const solution_t& up, const solution_t& down,
                                const solution_t& middle,
                                double y, double t) -> double {
                    return x_conditions->left_function({y}, t);
                };
            } else {
                functions[0] = [this, h_x, h_y]( const solution_t& left, const solution_t& right,
                                      const solution_t& up, const solution_t& down,
                                      const solution_t& middle,
                                      double y, double t) -> double {
                    const auto result =
                                   1.0 / c_function({ 0.0, y})
                                   * (1.0 / (h_x * h_x)
                                      * (lambda_function({ 0.5, y})
                                         * (right - middle)
                                         - lambda_function({ -0.5, y})
                                           * ( middle
                                           - ( right - 2.0 * h_x * x_conditions->left_function( {y}, t))))
                                      + 1.0 / (h_y * h_y)
                                        * (lambda_function({ 0.0, y + 0.5})
                                           * ( up - middle)
                                           - lambda_function({ 0.0, y - 0.5})
                                             * (middle - down)))
                                   + f_function({ 0, y}, t);
                    return result;
                };
            }
            if (!x_conditions->is_right_derivative()) {
                functions[1] = [this]( const solution_t& left, const solution_t& right,
                                 const solution_t& up, const solution_t& down,
                                 const solution_t& middle,
                                 double y, double t) -> double {
                    return x_conditions->right_function({y}, t);
                };
            } else {
                functions[1] = [this, h_x, h_y]( const solution_t& left, const solution_t& right,
                                           const solution_t& up, const solution_t& down,
                                           const solution_t& middle,
                                           double y, double t) -> double {
                    const auto result =
                            1.0 / c_function({ 0.0, y})
                            * (1.0 / (h_x * h_x)
                               * (lambda_function({ 0.5, y})
                                  * ( ( left + 2.0 * h_x * x_conditions->right_function({y},t)) - middle)
                                  - lambda_function({ -0.5, y})
                                    * ( middle - left))
                               + 1.0 / (h_y * h_y)
                                 * (lambda_function({ 0.0, y + 0.5})
                                    * ( up - middle)
                                    - lambda_function({ 0.0, y - 0.5})
                                      * (middle - down)))
                            + f_function({ 0, y}, t);
                    return result;
                };
            }
            if (!y_conditions->is_left_derivative()) {
                functions[2] = [this](const solution_t& left, const solution_t& right,
                                const solution_t& up, const solution_t& down,
                                const solution_t& middle,
                                double x, double t) -> double {
                    return y_conditions->left_function({x}, t);
                };
            } else {
                functions[2] = [this, h_x, h_y](const solution_t& left, const solution_t& right,
                                          const solution_t& up, const solution_t& down,
                                          const solution_t& middle,
                                          double x, double t) -> double {
                    const auto result =
                            1.0 / c_function({0.0, x})
                            * (1.0 / (h_x * h_x)
                               * (lambda_function({0.5, x})
                                  * (right - middle)
                                  - lambda_function({-0.5, x})
                                    * ( middle - right))
                               + 1.0 / (h_y * h_y)
                                 * (lambda_function({0.0, x + 0.5})
                                    * ( up - middle)
                                    - lambda_function({0.0, x - 0.5})
                                      * (middle - ( up - 2.0 * h_y * y_conditions->left_function({x}, t)))))
                            + f_function({0, x}, t);
                    return result;
                };
            }
            if (!y_conditions->is_left_derivative()) {
                functions[3] = [this](const solution_t& left, const solution_t& right,
                                 const solution_t& up, const solution_t& down,
                                 const solution_t& middle,
                                 double x, double t) -> double {
                    return y_conditions->right_function({x}, t);
                };
            } else {
                functions[3] = [this, h_x, h_y]( const solution_t& left, const solution_t& right,
                                           const solution_t& up, const solution_t& down,
                                           const solution_t& middle,
                                           double y, double t) -> double {
                    const auto result =
                            1.0 / c_function({ 0.0, y})
                            * (1.0 / (h_x * h_x)
                               * (lambda_function({ 0.5, y})
                                  * (right - middle)
                                  - lambda_function({ -0.5, y})
                                    * ( middle - left ))
                               + 1.0 / (h_y * h_y)
                                 * (lambda_function({ 0.0, y + 0.5})
                                    * ( down + 2.0 * y_conditions->right_function({y},t) - middle)
                                    - lambda_function({ 0.0, y - 0.5})
                                      * (middle - down)))
                            + f_function({ 0, y}, t);
                    return result;
                };
            }

            const auto left_x = functions[0];
            const auto right_x = functions[1];
            const auto left_y = functions[2];
            const auto right_y = functions[3];
            const auto ode_function = [this, h_x, h_y, left_x, left_y, right_x, right_y]
                    (const matrix<heat_equation_solver::solution_t> &old_solution, double t)
                    -> matrix<heat_equation_solver::solution_t> {
                auto result = old_solution;
                result[0][0] = left_x( result[0][0],result[0][0],result[0][0],result[0][0],result[0][0], 0, t );
                result[0][result.width() - 1] = left_x( result[0][0],result[0][0],result[0][0],result[0][0],result[0][0], (result.width() - 1) * h_y, t );
                result[result.height() - 1][0] = left_x( result[0][0],result[0][0],result[0][0],result[0][0],result[0][0], 0, t );
                result[result.height() - 1][result.width() - 1] = left_x( result[0][0],result[0][0],result[0][0],result[0][0],result[0][0], (result.width() - 1) * h_y, t );
                for (int i = 1; i < result.width() - 1; ++i) {
                    result[0][i] = left_x( result[0][i],
                                           result[1][i],
                                           result[0][i + 1],
                                           result[0][i - 1],
                                           result[0][i],
                                           i * h_y,
                                           t);
                    const auto N = result.height() - 1;
                    result[N][i] = right_x(result[N][i],
                                           result[N][i],
                                           result[N][i + 1],
                                           result[N][i - 1],
                                           result[N][i],
                                           i * h_y,
                                           t);
                }
                for (int i = 1; i < result.height() - 1; ++i) {
                    result[i][0] = left_y( result[i - 1][0],
                                           result[i + 1][0],
                                           result[i][1],
                                           result[i][0],
                                           result[i][0],
                                           i * h_x,
                                           t);
                    const auto N = result.width() - 1;
                    result[i][N] = right_y(result[i - 1][N],
                                           result[i + 1][N],
                                           result[i][N],
                                           result[i][N - 1],
                                           result[i][N],
                                           i * h_x,
                                           t);
                }

                for (int i = 1; i < old_solution.height() - 1; ++i) {
                    for (int j = 1; j < old_solution.width() - 1; ++j) {
                        const auto middle_point = old_solution[i][j];
                        result[i][j] =
                                1.0 / c_function({i * h_x, j * h_y})
                                * (1.0 / (h_x * h_x)
                                   * (lambda_function({i * h_x + 0.5, j * h_y})
                                      * (old_solution[i + 1][j] - middle_point)
                                      - lambda_function({i * h_x - 0.5, j * h_y})
                                        * (middle_point - old_solution[i - 1][j]))
                                   + 1.0 / (h_y * h_y)
                                     * (lambda_function({i * h_x, j * h_y + 0.5})
                                        * (old_solution[i][j + 1] - middle_point)
                                        - lambda_function({i * h_x, j * h_y - 0.5})
                                          * (middle_point - old_solution[i][j - 1])))
                                + f_function({i * h_x, j * h_y}, t);
                    }
                }
                return result;
            };
            ode_solver.set_function( ode_function );
        };


        return_t<matrix<solution_t>> get_current_layer() {
            const auto result = ode_solver.get_last_element();
            return {result, 0.0};
        }

        return_t<matrix<solution_t>> get_next_step() {
            const auto result = ode_solver.get_next_step( time_step );
            time += time_step;
            return result;
        }
    };

}
#endif //MAIN_CPP_HEAT_EQUATION_2D_SOLVER_HPP
