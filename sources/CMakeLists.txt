cmake_minimum_required(VERSION 3.16)
project(project_MMOOMPH)

set(CMAKE_CXX_STANDARD 20)

add_executable(folder main.cpp solution_t_concept.hpp
        vector.hpp  vector.cpp
		matrix.hpp
		SNAE_solver.hpp multisteps_methods.hpp one_step_methods.hpp ring_buffer.hpp
        tests.hpp
        SNAE_solver_tests.cpp SNAE_solver_tests.hpp
        matrix_tests.cpp matrix_tests.hpp vector_tests.cpp vector_tests.hpp 
		one_step_methods_tests.cpp one_step_methods_tests.hpp
		min_and_max_functions.hpp min_and_max_functions.cpp
		multisteps_methods_tests.cpp multisteps_methods_tests.hpp
		ring_buffer_tests.cpp ring_buffer_tests.hpp
		ODE_solver.hpp
		boundary_value_problem.hpp heat_equation_solver.hpp)