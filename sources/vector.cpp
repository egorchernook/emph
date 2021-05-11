#include "vector.hpp"

template<>
Numerical_methods::vector<double> Numerical_methods::inverse<double>( const Numerical_methods::vector<double> &initial){
    Numerical_methods::vector<double> result = initial;
    for( auto &x : result){
        x = 1.0/x;
    }
    return result;
}

double Numerical_methods::inverse( double initial){
    return 1.0 / initial;
}