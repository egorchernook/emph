#include "matrix.hpp"

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

    matrix<double> inverse_Gauss(const matrix<double> &matrix_) {

        assert( matrix_.height() == matrix_.width()) ;
        const int Height = matrix_.height();

        matrix<double> identity_matrix(0.0, Height);

        for (int i = 0; i < Height; i++)
            for (int j = 0; j < Height; j++)
                identity_matrix[i][j] = i == j ? 1.0 : 0.0;

        matrix<double> temporary_matrix = matrix_;

        for (int k = 0; k < Height; k++) {
            for (int i = k; i < Height; i++) {
                const auto multiplier = std::abs( temporary_matrix[i][k]) < std::numeric_limits<double>::epsilon() * 10'000
                                        ? 1.0 : temporary_matrix[i][k];
                for (int j = 0; j < Height; j++) {
                    temporary_matrix[i][j] *= inverse(multiplier);
                    identity_matrix[i][j] *= inverse(multiplier);
                }
            }

            for (int i = k + 1; i < Height; i++) {
                for (int j = 0; j < Height; j++) {
                    temporary_matrix[i][j] -= temporary_matrix[k][j];
                    identity_matrix[i][j] -= identity_matrix[k][j];
                }
            }
        }

        for (int i = 0; i < Height; i++)
            assert(temporary_matrix[i][i] != 0);

        for (int k = Height - 1; k >= 0; k--) {
            for (int i = k - 1; i >= 0; i--) {
                const auto multiplier = temporary_matrix[i][k];
                for (int j = Height - 1; j >= 0; j--) {
                    temporary_matrix[i][j] -= temporary_matrix[k][j] * multiplier;
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
            const matrix<double> &matrix_,
            const vector<double> &free_vector) {

        assert( matrix_.height() == matrix_.width());
        assert( matrix_.height() == free_vector.size());
        const int Height = matrix_.height();

        vector<double> temporary_free_vector = free_vector;
        matrix<double> temporary_matrix = matrix_;

        for (int k = 0; k < Height; k++) {
            for (int i = k; i < Height; i++) {
                const auto multiplier = std::abs( temporary_matrix[i][k]) < std::numeric_limits<double>::epsilon() * 10'000
                                        ? 1.0 : temporary_matrix[i][k];
                for (int j = 0; j < Height; j++) {
                    temporary_matrix[i][k] *= inverse(multiplier);
                }
                temporary_free_vector[i] *= inverse(multiplier);
            }

            for (int i = k + 1; i < Height; i++) {
                for (int j = 0; j < Height; j++) {
                    temporary_matrix[i][j] -= temporary_matrix[k][j];
                }
                temporary_free_vector[i] -= temporary_free_vector[k];
            }
        }

        vector<double> x(temporary_free_vector.size(), 0.0);
        for (int k = Height - 1; k >= 0; k--) {
            double sum{0.0};
            for (int i = k + 1; i < Height; i++) {
                sum += x[i] * temporary_matrix[k][i];
            }
            x[k] = (temporary_free_vector[k] - sum) * inverse(temporary_matrix[k][k]);

        }
        return x;
    }
}