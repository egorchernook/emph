#ifndef MAIN_CPP_MATRIX_HPP
#define MAIN_CPP_MATRIX_HPP

#include "vector.hpp"
#include "solution_t_concept.hpp"

namespace Numerical_methods {

    template<solution_element_t value_t, std::size_t Height, std::size_t Width = Height>
    class matrix {
    private:
        value_t ** data = nullptr;

    public:
        explicit matrix(value_t value) {

            data = new value_t *[Height];

            for (int i = 0; i < Height; ++i) {
                data[i] = new value_t[Width];
                for (int j = 0; j < Width; j++) {
                    data[i][j] = value;
                }
            }
        }

        matrix(){
            data = new value_t *[Height];
            for (int i = 0; i < Height; ++i) {
                data[i] = new(std::nothrow) value_t[Width];
            }
        }

        ~matrix() noexcept {
            if ( data != nullptr ) {
                for (int i = 0; i < Height; ++i) {
                    if (data[i] != nullptr) {
                        delete[] data[i];
                    }
                }
                delete[] data;
                data = nullptr;
            }
        }

        matrix(
                const std::initializer_list<
                        std::initializer_list<
                                value_t>> &list)
        {
            assert( Height == list.size());
            for (auto &x : list) {
                assert(Width == x.size());
            }
            data = new value_t *[Height];
            auto it = list.begin();
            for (int i = 0; i < Height; ++i, ++it) {
                data[i] = new value_t[Width];
                std::move(it->begin(), it->end(), data[i]);
            }
        }

        matrix( const matrix& other){
            data = new value_t *[Height];
            for ( int i = 0; i < Height; i++){
                data[i] = new(std::nothrow) value_t[Width];
                for (int j = 0; j < Width; j++)
                    data[i][j] = other.data[i][j];
            }
        }

        [[nodiscard]] constexpr std::size_t height() const noexcept {
            return Height;
        }

        [[nodiscard]] constexpr std::size_t width() const noexcept {
            return Width;
        }

        [[nodiscard]] value_t * operator[](std::size_t idx) const {
            assert( idx < Height );
            return data[idx];
        }

        [[nodiscard]] value_t * at(std::size_t idx) const {
            assert( idx < Height );
            return data[idx];
        }

        [[nodiscard]] value_t& at(std::size_t idx, std::size_t idy) const {
            assert(idx < Height);
            assert(idy < Width);
            return data[idx][idy];
        }

        matrix& operator=(const matrix &other) {
            if (this == &other)
                return *this;

            for (int i = 0; i < Height; ++i) {
                for (int j = 0; j < Width; ++j) {
                    data[i][j] = other[i][j];
                }
            }
            return *this;
        }

        matrix<value_t, Height, Width>& operator+=(const matrix <value_t, Height, Width> &another) {

            assert(this->height() == another.height());
            assert(this->width() == another.width());

            for (int idx_x = 0; idx_x < Height; idx_x++) {
                for (int idx_y = 0; idx_y < Width; idx_y++) {
                    this->at(idx_x, idx_y) += another[idx_x][idx_y];
                }
            }
            return *this;
        }

        matrix<value_t, Height, Width> &operator-=(const matrix <value_t, Height, Width> &another) {

            assert(this->height() == another.height());
            assert(this->width() == another.width());

            for (int idx_x = 0; idx_x < Height; idx_x++) {
                for (int idx_y = 0; idx_y < Width; idx_y++) {
                    this->at(idx_x, idx_y) -= another[idx_x][idx_y];
                }
            }
            return *this;
        }

        matrix<value_t, Height, Width> &operator*=(double value) {

            for (int idx_x = 0; idx_x < Height; idx_x++) {
                for (int idx_y = 0; idx_y < Width; idx_y++) {
                    this->at(idx_x, idx_y) *= value;
                }
            }
            return *this;
        }

        matrix<value_t, Height, Width> &operator/=(double value) {

            for (int idx_x = 0; idx_x < Height; idx_x++) {
                for (int idx_y = 0; idx_y < Width; idx_y++) {
                    this->at(idx_x, idx_y) /= value;
                }
            }
            return *this;
        }

    };

    template<solution_element_t value_t, std::size_t Height, std::size_t Width = Height>
    matrix<value_t, Height, Width> operator+(
            const matrix <value_t, Height, Width> &lhs,
            const matrix <value_t, Height, Width> &rhs) {

        matrix<value_t, Height, Width> result;
        for (int idx_x = 0; idx_x < Height; idx_x++) {
            for (int idx_y = 0; idx_y < Width; idx_y++) {
                result[idx_x][idx_y] = lhs[idx_x][idx_y] + rhs[idx_x][idx_y];
            }
        }
        return result;
    }

    template<solution_element_t value_t, std::size_t Height, std::size_t Width = Height>
    matrix<value_t, Height, Width> operator-(
            const matrix <value_t, Height, Width> &lhs,
            const matrix <value_t, Height, Width> &rhs) {

        matrix<value_t, Height, Width> result;
        for (int idx_x = 0; idx_x < Height; idx_x++) {
            for (int idx_y = 0; idx_y < Width; idx_y++) {
                result[idx_x][idx_y] = lhs[idx_x][idx_y] - rhs[idx_x][idx_y];
            }
        }
        return result;
    }

    template<solution_element_t value_t, std::size_t Height, std::size_t Width = Height>
    matrix<value_t, Height, Width> operator*(
            matrix<value_t, Height, Width> &lhs,
            double rhs) {

        matrix<value_t, Height, Width> result;
        for (int idx_x = 0; idx_x < Height; idx_x++) {
            for (int idx_y = 0; idx_y < Width; idx_y++) {
                result[idx_x][idx_y] = lhs[idx_x][idx_y] * rhs;
            }
        }
        return result;
    }

    template<solution_element_t value_t, std::size_t Height, std::size_t Width = Height>
    matrix<value_t, Height, Width> operator*(
            double lhs,
            const matrix <value_t, Height, Width> &rhs) {

        matrix<value_t, Height, Width> result;
        for (int idx_x = 0; idx_x < Height; ++idx_x) {
            for (int idx_y = 0; idx_y < Width; ++idx_y) {
                result[idx_x][idx_y] = lhs * rhs[idx_x][idx_y];
            }
        }
        return result;
    }

    template<solution_element_t value_t, std::size_t lhs_Height, std::size_t lhs_Width, std::size_t rhs_Height, std::size_t rhs_Width>
    matrix<value_t, lhs_Height, rhs_Width> operator*(
            const matrix<value_t, lhs_Height, lhs_Width>& lhs,
            const matrix<value_t, rhs_Height, rhs_Width>& rhs)
    {
        assert( lhs_Width == rhs_Height );

        matrix<value_t,lhs_Height, rhs_Width> buf{};

        for (int i = 0; i < lhs_Height; ++i)
        {
            for (int j = 0; j < rhs_Width; ++j)
            {
                value_t temp = lhs[i][0] * rhs[0][j];
                for (int k = 1; k < lhs_Width; ++k)
                {
                    temp += lhs[i][k] * rhs[k][j];
                }
                buf[i][j] = temp;
            }
        }
        return buf;
    }

    template<solution_element_t value_t, std::size_t Height, std::size_t Width = Height>
    matrix<value_t, Height, Width> operator/(
            const matrix <value_t, Height, Width> &lhs,
            double rhs) {

        matrix<value_t, Height, Width> result;
        for (int idx_x = 0; idx_x < Height; idx_x++) {
            for (int idx_y = 0; idx_y < Width; idx_y++) {
                result[idx_x][idx_y] = lhs[idx_x][idx_y] / rhs;
            }
        }
        return result;
    }

    template<std::size_t Height, std::size_t Width = Height>
    std::ostream &operator<<(std::ostream &stream, const matrix<double, Height, Width> &output_matrix) {
        const int length = 10;
        for (int i = 0; i < Height; i++) {
            for (int j = 0; j < Width; j++) {
                stream.width(length);
                stream << output_matrix[i][j] << "\t";
            }
            stream << std::endl;
        }
        return stream;
    }

    template<std::size_t Height>
    matrix<double, Height, Height> inverse_Gauss(const matrix<double, Height, Height> &matrix_) {

        matrix<double, Height> identity_matrix{double()};
        for (int i = 0; i < Height; i++)
            for (int j = 0; j < Height; j++)
                identity_matrix[i][j] = i == j ? double{1.0} : double{0.0};

        matrix<double, Height> temporary_matrix = matrix_;
        double multiplier = 0.0;

        for (int k = 0; k < Height; k++) {
            for (int i = k; i < Height; i++) {
                multiplier = temporary_matrix[i][k] == 0.0 ? 1.0 : temporary_matrix[i][k];
                for (int j = 0; j < Height; j++) {
                    temporary_matrix[i][j] /= multiplier;
                    identity_matrix[i][j] /= multiplier;
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
                multiplier = temporary_matrix[i][k];
                for (int j = Height - 1; j >= 0; j--) {
                    temporary_matrix[i][j] -= temporary_matrix[k][j] * multiplier;
                    identity_matrix[i][j] -= identity_matrix[k][j] * multiplier;
                }
            }
        }

        return identity_matrix;
    }

    template<std::size_t Height>
    vector<double> Gauss_method(
            const matrix<double, Height, Height> &matrix_,
            const vector<double> &free_vector) {

        vector<double> temporary_free_vector = free_vector;
        matrix<double, Height> temporary_matrix = matrix_;
        double multiplier{0.0};

        for (int k = 0; k < Height; k++) {
            for (int i = k; i < Height; i++) {
                multiplier = temporary_matrix[i][k] == 0.0 ? 1.0 : temporary_matrix[i][k];
                for (int j = 0; j < Height; j++) {
                    temporary_matrix[i][k] /= multiplier;
                }
                temporary_free_vector[i] /= multiplier;
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
            double sum = 0.0;
            for (int i = k + 1; i < Height; i++) {
                sum += x[i] * temporary_matrix[k][i];
            }
            x[k] = (temporary_free_vector[k] - sum) / temporary_matrix[k][k];

        }
        return x;
    }

}

#endif //MAIN_CPP_MATRIX_HPP
