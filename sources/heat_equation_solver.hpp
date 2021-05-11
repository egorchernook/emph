#ifndef MAIN_CPP_HEAT_EQUATION_SOLVER_HPP
#define MAIN_CPP_HEAT_EQUATION_SOLVER_HPP

#include <functional>
#include <utility>
#include "solution_t_concept.hpp"

namespace Numerical_methods{

    template<solution_element_t solution_t, std::size_t dimension>
    struct boundary_conditions {
        using coordinate_t = vector<double>;
        using function_t = std::function< coordinate_t(const coordinate_t&, double)>;
        const function_t left_value;
        const function_t right_value;
        const coordinate_t left_border;
        const coordinate_t right_border;
        boundary_conditions(
                            function_t left_value_,
                            function_t right_value_,
                            coordinate_t left_border_,
                            coordinate_t right_border_)
                            :
                            left_value(left_value_),
                            right_value(right_value_),
                            left_border(left_border_),
                            right_border(right_border_)
                            {
            assert( left_border.size() == dimension - 1);
            assert( right_border.size() == dimension - 1);
        }
    };

    template<solution_element_t solution_t, std::size_t size>
    class Dirichlet_conditions : public boundary_conditions<solution_t, size>{
        using boundary_conditions<solution_t, size>::boundary_conditions();
    };

    template<solution_element_t solution_t, std::size_t size>
    class Neumann_conditions : public boundary_conditions<solution_t, size>{
        using boundary_conditions<solution_t, size>::boundary_conditions();
    };

    template<solution_element_t solution_t, std::size_t size>
    class mixed_conditions : public boundary_conditions<solution_t, size>{
        using boundary_conditions<solution_t, size>::boundary_conditions();
    };

    template<solution_element_t solution_t, std::size_t size>
    class first_mixed_conditions : public mixed_conditions<solution_t, size>{
        using mixed_conditions<solution_t, size>::mixed_conditions();
    };

    template<solution_element_t solution_t, std::size_t size>
    class second_mixed_conditions : public mixed_conditions<solution_t, size>{
        using mixed_conditions<solution_t, size>::mixed_conditions();
    };

    template<std::size_t dimension>
    class heat_equation_solver {
    private:
        std::function< double( vector<double>)> initial_conditions;
    public:
        const boundary_conditions<double,dimension> *const conditions = nullptr;
        heat_equation_solver( std::function< double( vector<double>)> initial_conditions_,
                              boundary_conditions<double,dimension> conditions_)
                              : initial_conditions(initial_conditions_){
            conditions = conditions_;
        }
        heat_equation_solver() noexcept {
            conditions = nullptr;
        }
    };

}
#endif //MAIN_CPP_HEAT_EQUATION_SOLVER_HPP
