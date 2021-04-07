#ifndef MAIN_CPP_RING_BUFFER_TESTS_HPP
#define MAIN_CPP_RING_BUFFER_TESTS_HPP

#include <vector>
#include "solution_t_concept.hpp"
#include "ring_buffer.hpp"

namespace Numerical_methods {
    class ring_buffer_tests {
    private:
    public:
        ring_buffer_tests();

        static bool push_test();

        static bool read_test();

        static bool current_state_test();
    };
}

#endif //MAIN_CPP_RING_BUFFER_TESTS_HPP
