#include <iostream>
#include "SNAE_solver_tests.hpp"

Numerical_methods::SNAE_solver_tests::SNAE_solver_tests() {
    std::cout << std::endl << "\ttesting solving methods for systems of nonlinear algebraic equations..." << std::endl;
    assert( test_fixed_point_iterations_method() );
    std::cout << "\t\tfixed point iterations method works" << std::endl;
    assert( test_Seidels_method() );
    std::cout << "\t\tSeidels method works" << std::endl;
    //assert( test_Newtons_method() );
    //std::cout << "\t\tNewtons method works" << std::endl;
}

bool Numerical_methods::SNAE_solver_tests::test_fixed_point_iterations_method() {
    vector<double> initial_state = {0.0,0.0};
    fixed_point_iterations_method< vector<double>> solver(SNAE_solver_tests::test_equation);
    auto solving_result = solver.solve( initial_state);

    bool result = ( solving_result.solution[0] - solving_result.error ) * 10'000 < 0.5 * 10'000 &&
                  ( solving_result.solution[0] + solving_result.error ) * 10'000 > 0.5 * 10'000 &&
                  ( solving_result.solution[1] - solving_result.error ) * 10'000 < 0.5 * 10'000 &&
                  ( solving_result.solution[1] + solving_result.error ) * 10'000 > 0.5 * 10'000;
    return result;
}

bool Numerical_methods::SNAE_solver_tests::test_Seidels_method() {
    vector<double> initial_state = {0.0, 0.0};
    Seidels_method<vector<double>> solver(SNAE_solver_tests::test_equation);
    auto solving_result = solver.solve( initial_state);

    bool result = (solving_result.solution[0] - solving_result.error) * 10'000 < 0.5 * 10'000 &&
                  (solving_result.solution[0] + solving_result.error) * 10'000 > 0.5 * 10'000 &&
                  (solving_result.solution[1] - solving_result.error) * 10'000 < 0.5 * 10'000 &&
                  (solving_result.solution[1] + solving_result.error) * 10'000 > 0.5 * 10'000;
    return result;
}

bool Numerical_methods::SNAE_solver_tests::test_Newtons_method() {
    vector<double> initial_state = {1.0, 0.0};
    Newtons_method<vector<double>> solver([](const vector<double>& old) -> vector<double>{
                                                                    return old - test_equation(old);
                                                                },
                                          { SNAE_solver_tests::test_equation_derivative_1,
                                                             SNAE_solver_tests::test_equation_derivative_2});
    auto solving_result = solver.solve( initial_state);

    bool result = (solving_result.solution[0] - solving_result.error) * 10'000 < 0.5 * 10'000 &&
                  (solving_result.solution[0] + solving_result.error) * 10'000 > 0.5 * 10'000 &&
                  (solving_result.solution[1] - solving_result.error) * 10'000 < 0.5 * 10'000 &&
                  (solving_result.solution[1] + solving_result.error) * 10'000 > 0.5 * 10'000;
    return result;
}
