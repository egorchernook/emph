#include <iostream>
#include <cassert>
#include <functional>
#include <algorithm>
#include <vector>

#include "solution_t_concept.hpp"
#include "ring_buffer.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "one_step_methods.hpp"
#include "multisteps_methods.hpp"

#include "tests.hpp"

constexpr int N = 4;
constexpr int omegas_amount = 10'000 ;
constexpr int steps_amount = 10'000;

int main(){
/*
    const Numerical_methods::implicit_Adams_Moulton_method< precision_order > impl_msm_solver{ {3.0/8, 19.0/24, -5.0/24, 1.0/24} };
    const Numerical_methods::implicit_Runge_Kutta_method< precision_order - 1> impl_rkm_solver{
            { { 1.0/6 , 0.0   , -1.0/6 },
              { 1.0/12, 5.0/12,  0.0   },
              { 0.5   , 1.0/3 ,  1.0/6 } },
            {   1.0/6 , 2.0/3 , 1.0/6    },
            {   0.0   , 1.0/2   , 1.0      }
    };
*/
    Numerical_methods::test();
    return 0;
}