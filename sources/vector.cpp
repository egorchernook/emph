#include "vector.hpp"

Numerical_methods::vector<double> Numerical_methods::inverse(const Numerical_methods::vector<double> &initial){
    Numerical_methods::vector<double> result = initial;
    for( auto &x : result){
        x = 1.0 / x;
    }
    return result;
}

double Numerical_methods::inverse( double initial){
    if ( std::abs( initial) < std::numeric_limits<double>::epsilon() * 10'000){
        return 1.0 / std::numeric_limits<double>::epsilon() * 10'000;
    } else {
        return 1.0 / initial;
    }
}