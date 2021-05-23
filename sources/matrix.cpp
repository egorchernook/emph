#include "matrix.hpp"
#include "min_and_max_functions.hpp"

namespace Numerical_methods{

    std::ostream &operator<<(std::ostream &stream, const matrix<double> &output_matrix) {
        const int length = 10;
        for (int i = 0; i < output_matrix.height(); i++) {
            for (int j = 0; j < output_matrix.width(); j++) {
                stream.width(length);
                stream << output_matrix[i][j] << "\t";
            }
            stream << std::endl;
        }
        return stream;
    }

    matrix<double> inverse_Gauss( matrix<double> matrix_) {

        assert( matrix_.height() == matrix_.width()) ;
        const int Height = matrix_.height();

        matrix<double> identity_matrix(0.0, Height);

        for (int i = 0; i < Height; i++)
            for (int j = 0; j < Height; j++)
                identity_matrix[i][j] = i == j ? 1.0 : 0.0;

        for (int k = 0; k < Height; k++) {
            for (int i = k; i < Height; i++) {
                const auto multiplier = std::abs( matrix_[i][k]) < std::numeric_limits<double>::epsilon() * 10'000
                                        ? 1.0 : matrix_[i][k];
                for (int j = 0; j < Height; j++) {
                    matrix_[i][j] /= multiplier;
                    identity_matrix[i][j] /= multiplier;
                }
            }

            for (int i = k + 1; i < Height; i++) {
                for (int j = 0; j < Height; j++) {
                    matrix_[i][j] -= matrix_[k][j];
                    identity_matrix[i][j] -= identity_matrix[k][j];
                }
            }
        }

        for (int i = 0; i < Height; i++)
            assert(matrix_[i][i] != 0);

        for (int k = Height - 1; k >= 0; k--) {
            for (int i = k - 1; i >= 0; i--) {
                const auto multiplier = matrix_[i][k];
                for (int j = Height - 1; j >= 0; j--) {
                    matrix_[i][j] -= matrix_[k][j] * multiplier;
                    identity_matrix[i][j] -= identity_matrix[k][j] * multiplier;
                }
            }
        }

        return identity_matrix;
    }

    matrix<double> inverse(const matrix<double> &matrix_){
        return inverse_Gauss(matrix_);
    }

    vector<double> Gauss_method(
            matrix<double> matrix_,
            vector<double> free_vector) {

        assert( matrix_.height() == matrix_.width());
        assert( matrix_.height() == free_vector.size());
        const int Height = static_cast<int>(matrix_.height());

        vector<int> result_indexes(free_vector.size());
        for( int i = 0; i < free_vector.size(); ++i){
            result_indexes[i] = i;
        }

        for (int k = 0; k < Height; k++) {

            double max_element = matrix_[k][k];
            int i_max{k};
            int j_max{k};
            for( int i = k; i < Height; ++i){
                for( int j = k; j < Height; ++j){
                    if( matrix_[i][j] > max_element) {
                        max_element = matrix_[i][j];
                        j_max = j;
                        i_max = i;
                    }
                }
            }

            if( i_max != k ) {
                const auto temp_ptr = matrix_[i_max];
                matrix_[i_max] = matrix_[k];
                matrix_[k] = temp_ptr;
                const auto temp_d = free_vector[i_max];
                free_vector[i_max] = free_vector[k];
                free_vector[k] = temp_d;

                if( j_max != k ) {
                    for( int i = 0; i < Height; ++i ) {
                        const double temp = matrix_[j_max][i];
                        matrix_[j_max][i] = matrix_[k][i];
                        matrix_[k][i] = temp;
                    }
                    result_indexes[k] = j_max;
                    result_indexes[j_max] = j_max;
                }
            }

            for (int i = k; i < Height; i++) {
                const auto multiplier =
                        std::abs(
                                matrix_[i][k]) <
                                std::numeric_limits<double>::epsilon() * 10'000
                                    ? 1.0 : matrix_[i][k];
                for (int j = 0; j < Height; j++) {
                    matrix_[i][k] /= multiplier;
                }
                free_vector[i] /= multiplier;
            }

            for (int i = k + 1; i < Height; i++) {
                for (int j = 0; j < Height; j++) {
                    matrix_[i][j] -= matrix_[k][j];
                }
                free_vector[i] -= free_vector[k];
            }
        }

        vector<double> x(free_vector.size(), 0.0);
        for (int k = Height - 1; k >= 0; k--) {
            double sum{0.0};
            for (int i = k + 1; i < Height; i++) {
                const auto temp = (
                                    std::abs(matrix_[k][i]) < std::numeric_limits<double>::epsilon() * 10'000 )
                                            ? std::numeric_limits<double>::epsilon() * 10'000 : matrix_[k][i];
                sum += x[i] * temp;
            }
            const auto temp = (
                                      std::abs(matrix_[k][k]) < std::numeric_limits<double>::epsilon() * 10'000 )
                              ? std::numeric_limits<double>::epsilon() * 10'000 : matrix_[k][k];
            x[k] = (free_vector[k] - sum) / temp;

        }

        vector<double> result = x;
        for( int i = 0; i < x.size(); ++i){
            result[ result_indexes[i] ] = x[i];
        }
        return result;
    }
}