#include "expand_function.hpp"

namespace Numerical_methods{

    matrix<double> expand( const matrix<double>& matrix){
        return matrix;
    }

    vector<double> expand( const vector<double>& vector){
        return vector;
    }

    template<>
    matrix<double> expand<vector<double>>( const vector<vector<double>>& vector ) {
        matrix<double> result(vector.size(), vector[0].size());
        for (int i = 0; i < result.height(); ++i) {
            for (int j = 0; j < result.width(); ++j) {
                result[i][j] = vector[i][j];
            }
        }
        return result;
    }
}