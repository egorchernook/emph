#include "vector.hpp"

namespace Numerical_methods {
    template<>
    vector<double> inverse<double>( const vector<double> &initial){
        vector<double> result = initial;
        for( auto &x : result){
            x = 1.0/x;
        }
        return result;
    }
}