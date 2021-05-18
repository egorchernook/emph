#ifndef MAIN_CPP_EXPAND_FUNCTION_HPP
#define MAIN_CPP_EXPAND_FUNCTION_HPP

#include "solution_t_concept.hpp"
#include "vector.hpp"
#include "matrix.hpp"

namespace Numerical_methods{

    matrix<double> expand( const matrix<double>&);

    vector<double> expand( const vector<double>&);

    template<solution_t_concept solution_t>
    matrix<double> expand( const matrix<solution_t>& matrix_){
        const auto temp = expand(matrix_[0][0]);
        matrix<double> result( matrix_.height() * temp.height(), matrix_.width() * temp.width());
        for( int i = 0; i < matrix_.height(); ++i ) {
            for( int j = 0; j < matrix_.width(); ++j ) {
                const auto inner = expand(matrix_[i][j]);
                for(int n = 0; n < inner.height(); ++n ) {
                    for(int m = 0; m < inner.width(); ++m ) {
                        result[ i * matrix_.height() + n][ j * matrix_.width() + m] = inner[n][m];
                    }
                }
            }
        }
        return result;
    }

    template<solution_t_concept solution_t>
    matrix<double> expand( const vector<solution_t>& vector){
        const auto temp = expand(vector[0]);
        matrix<double> result( vector.size() * temp.height(), temp.width());
        for( int i = 0; i < vector.size(); ++i ) {
            const auto inner = expand(vector[i]);
            for (int n = 0; n < inner.height(); ++n) {
                for (int m = 0; m < inner.width(); ++m) {
                    result[i * vector.size() + n][m] = inner[n][m];
                }
            }
        }
        return result;
    }

    template<>
    matrix<double> expand<vector<double>>( const vector<vector<double>>& vector);
}

#endif //MAIN_CPP_EXPAND_FUNCTION_HPP
