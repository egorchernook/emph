#ifndef MAIN_CPP_MATRIX_HPP
#define MAIN_CPP_MATRIX_HPP

#include <vector>

#include "vector.hpp"
#include "solution_t_concept.hpp"

namespace Numerical_methods {

    template<typename value_t>
    class base_matrix{
    protected:
        value_t ** data = nullptr;
        std::size_t Height = 0;
        std::size_t Width = 0;
    public:
        base_matrix(value_t value, std::size_t height, std::size_t width = 0) : Height(height), Width(width){

            if( width == 0){
                Width = height;
            }
            data = new value_t *[Height];

            for (int i = 0; i < Height; ++i) {
                data[i] = new value_t[Width];
                for (int j = 0; j < Width; j++) {
                    data[i][j] = value;
                }
            }
        }

        base_matrix( std::size_t height, std::size_t width = 0) : Height(height), Width(width) {
            if( width == 0){
                Width = height;
            }
            data = new value_t *[Height];
            for (int i = 0; i < Height; ++i) {
                data[i] = new value_t[Width];
            }
        }

        ~base_matrix() noexcept {
            if (data != nullptr ) {
                for (int i = 0; i < Height; ++i) {
                    if (data[i] != nullptr) {
                        delete[] data[i];
                    }
                }
                delete[] data;
                data = nullptr;
            }
        }

        base_matrix(
                const std::initializer_list<
                        std::initializer_list<
                                value_t>> &list)
                : Height(list.size()), Width(list.begin()->size())
        {
            assert(Height == list.size());
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

        base_matrix( const base_matrix& other) : Height(other.height()), Width(other.width()) {
            data = new value_t *[Height];
            for (int i = 0; i < Height; i++){
                data[i] = new value_t[Width];
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

        [[nodiscard]] value_t *& operator[](std::size_t idx) const {
            assert(idx < Height );
            return data[idx];
        }

        [[nodiscard]] value_t *& at(std::size_t idx) const {
            assert(idx < Height );
            return data[idx];
        }

        [[nodiscard]] value_t& at(std::size_t idx, std::size_t idy) const {
            assert(idx < Height);
            assert(idy < Width);
            return data[idx][idy];
        }

        base_matrix& operator=(const base_matrix &other) {
            if (this == &other)
                return *this;
            delete [] data;
            Height = other.height();
            Width = other.width();
            data = new value_t *[Height];
            for (int i = 0; i < Height; i++){
                data[i] = new value_t[Width];
                for (int j = 0; j < Width; j++)
                    data[i][j] = other.data[i][j];
            }
            return *this;
        }
    };

    template<solution_element_t value_t>
    class matrix : public base_matrix<value_t>{
        using base_t = base_matrix<value_t>;
    public:
        using base_t::base_t;
        matrix<value_t>& operator+=(const matrix <value_t> &another) {
            assert(this->height() == another.height());
            assert(this->width() == another.width());

            for (int idx_x = 0; idx_x < base_t::Height; idx_x++) {
                for (int idx_y = 0; idx_y < base_t::Width; idx_y++) {
                    this->at(idx_x, idx_y) += another[idx_x][idx_y];
                }
            }
            return *this;
        }

        matrix<value_t> &operator-=(const matrix <value_t> &another) {
            assert(this->height() == another.height());
            assert(this->width() == another.width());

            for (int idx_x = 0; idx_x < base_t::Height; idx_x++) {
                for (int idx_y = 0; idx_y < base_t::Width; idx_y++) {
                    this->at(idx_x, idx_y) -= another[idx_x][idx_y];
                }
            }
            return *this;
        }

        matrix<value_t> &operator*=(double value) {
            for (int idx_x = 0; idx_x < base_t::Height; idx_x++) {
                for (int idx_y = 0; idx_y < base_t::Width; idx_y++) {
                    this->at(idx_x, idx_y) *= value;
                }
            }
            return *this;
        }

        matrix<value_t> &operator/=(double value) {
            for (int idx_x = 0; idx_x < base_t::Height; idx_x++) {
                for (int idx_y = 0; idx_y < base_t::Width; idx_y++) {
                    this->at(idx_x, idx_y) /= value;
                }
            }
            return *this;
        }

    };

    template<solution_element_t value_t>
    matrix<value_t> operator+(
            const matrix <value_t> &lhs,
            const matrix <value_t> &rhs) {

        assert(lhs.height() == rhs.height());
        assert(lhs.width() == rhs.width());

        matrix<value_t> result(lhs.height(), lhs.width());
        for (int idx_x = 0; idx_x < result.height(); ++idx_x) {
            for (int idx_y = 0; idx_y < result.width(); ++idx_y) {
                result[idx_x][idx_y] = lhs[idx_x][idx_y] + rhs[idx_x][idx_y];
            }
        }
        return result;
    }

    template<solution_element_t value_t>
    matrix<value_t> operator-(
            const matrix <value_t> &lhs,
            const matrix <value_t> &rhs) {

        assert(lhs.height() == rhs.height());
        assert(lhs.width() == rhs.width());
        matrix<value_t> result(lhs.height(), lhs.width());
        for (int idx_x = 0; idx_x < result.height(); ++idx_x) {
            for (int idx_y = 0; idx_y < result.width(); ++idx_y) {
                result[idx_x][idx_y] = lhs[idx_x][idx_y] - rhs[idx_x][idx_y];
            }
        }
        return result;
    }

    template<solution_element_t value_t>
    matrix<value_t> operator*(
            matrix<value_t> &lhs,
            double rhs) {

        matrix<value_t> result(lhs.height(), lhs.width());
        for (int idx_x = 0; idx_x < result.height(); ++idx_x) {
            for (int idx_y = 0; idx_y < result.width(); ++idx_y) {
                result[idx_x][idx_y] = lhs[idx_x][idx_y] * rhs;
            }
        }
        return result;
    }

    template<solution_element_t value_t>
    matrix<value_t> operator*(
            double lhs,
            const matrix <value_t> &rhs) {

        matrix<value_t> result(rhs.height(), rhs.width());
        for (int idx_x = 0; idx_x < result.height(); ++idx_x) {
            for (int idx_y = 0; idx_y < result.width(); ++idx_y) {
                result[idx_x][idx_y] = lhs * rhs[idx_x][idx_y];
            }
        }
        return result;
    }

    template<solution_element_t value_t>
    matrix<value_t> operator*(
            const matrix<value_t>& lhs,
            const matrix<value_t>& rhs)
    {
        assert( lhs.width() == rhs.height() );

        matrix<value_t> buf( lhs.width(), rhs.height());

        for (int i = 0; i < lhs.height(); ++i)
        {
            for (int j = 0; j < rhs.width(); ++j)
            {
                value_t temp = lhs[i][0] * rhs[0][j];
                for (int k = 1; k < lhs.width(); ++k)
                {
                    temp += lhs[i][k] * rhs[k][j];
                }
                buf[i][j] = temp;
            }
        }
        return buf;
    }

    template<solution_element_t value_t>
    matrix<value_t> operator/(
            const matrix <value_t> &lhs,
            double rhs) {

        matrix<value_t> result(lhs.height(), lhs.width());
        for (int idx_x = 0; idx_x < result.height(); ++idx_x) {
            for (int idx_y = 0; idx_y < result.width(); ++idx_y) {
                result[idx_x][idx_y] = lhs[idx_x][idx_y] / rhs;
            }
        }
        return result;
    }

    std::ostream &operator<<(std::ostream &stream, const matrix<double> &output_matrix);

    matrix<double> inverse_Gauss(matrix<double> matrix_);

    matrix<double> inverse(const matrix<double> &matrix_);

    vector<double> Gauss_method(
            matrix<double> matrix_,
            vector<double> free_vector);
}

#endif //MAIN_CPP_MATRIX_HPP
