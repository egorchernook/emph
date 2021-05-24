#ifndef MAIN_CPP_VECTOR_HPP
#define MAIN_CPP_VECTOR_HPP

#include <cstdlib>
#include <vector>

namespace Numerical_methods {
    #ifdef __has_include
        #if __has_include(<concepts>)
            #include <concepts>
            template<typename T>
            concept integral = std::integral<T>;
        #else
            #if defined __GNUC__ && __GNUC__ < 7
                #define integral typename
            #else
                template<typename T>
                concept integral = std::is_integral<T>::value;
            #endif
        #endif
    #endif
}

#include "solution_t_concept.hpp"

namespace Numerical_methods {

    template<solution_element_t T>
    class vector : private std::vector<T> {
    public:
        using std::vector<T>::vector;
        using std::vector<T>::operator[];
        using std::vector<T>::at;
        using std::vector<T>::push_back;
        using std::vector<T>::size;
        using std::vector<T>::operator=;
        using std::vector<T>::reserve;
        using std::vector<T>::begin;
        using std::vector<T>::end;
        using std::vector<T>::rbegin;
        using std::vector<T>::rend;
        using std::vector<T>::front;
        using std::vector<T>::back;

        template<integral value_t>
        explicit vector( std::size_t count, value_t value = 0.0) : std::vector<T>(count, T(value) ) {
            const T val{value};
            for( auto& x : *this){
                x = val;
            }
        }

        vector<T> &operator+=(const vector<T> &another) {

            assert(this->size() == another.size());

            for (std::size_t idx = 0; idx < this->size(); idx++) {
                this->at(idx) += another[idx];
            }
            return *this;
        }

        vector<T> &operator-=(const vector<T> &another) {

            assert(this->size() == another.size());

            for (std::size_t idx = 0; idx < this->size(); idx++) {
                this->at(idx) -= another[idx];
            }
            return *this;
        }

        vector<T> &operator*=(double value) {

            for (std::size_t idx = 0; idx < this->size(); idx++) {
                this->at(idx) *= value;
            }
            return *this;
        }

        vector<T> &operator/=(double value) {

            for (std::size_t idx = 0; idx < this->size(); idx++) {
                this->at(idx) /= value;
            }
            return *this;
        }

        vector<T> &operator*=(const vector<T>& another) {
            assert( this->size() == another.size() );
            for (std::size_t idx = 0; idx < this->size(); idx++) {
                this->at(idx) *= another[idx];
            }
            return *this;
        }
    };

    template<solution_element_t T>
    vector<T> operator+(
            const vector<T> &lhs,
            const vector<T> &rhs) {

        assert(lhs.size() == rhs.size());

        vector<T> result{};
        for (std::size_t idx = 0; idx < lhs.size(); idx++) {
            result.push_back(lhs[idx] + rhs[idx]);
        }
        return result;
    }

    template<solution_element_t T>
    vector<T> operator-(
            const vector<T> &lhs,
            const vector<T> &rhs) {

        assert(lhs.size() == rhs.size());

        vector<T> result{};
        for (std::size_t idx = 0; idx < lhs.size(); idx++) {
            result.push_back(lhs[idx] - rhs[idx]);
        }
        return result;
    }

    template<solution_element_t T>
    vector<T> operator*(
            const vector<T> &lhs,
            double rhs) {

        vector<T> result{};
        for (std::size_t idx = 0; idx < lhs.size(); idx++) {
            result.push_back(lhs[idx] * rhs);
        }
        return result;
    }

    template<solution_element_t T>
    vector<T> operator*(
            double lhs,
            const vector<T> &rhs) {

        vector<T> result{};
        for (std::size_t idx = 0; idx < rhs.size(); idx++) {
            result.push_back(lhs * rhs[idx]);
        }
        return result;
    }

    template<solution_element_t T>
    vector<T> operator/(
            const vector<T> &lhs,
            double rhs) {

        vector<T> result{};
        for (std::size_t idx = 0; idx < lhs.size(); idx++) {
            result.push_back(lhs[idx] / rhs);
        }
        return result;
    }

    template<solution_element_t T>
    vector<T> operator*(
            const vector<T> &lhs,
            const vector<T> &rhs) {

        assert( lhs.size() == rhs.size() );
        vector<T> result{};
        for (std::size_t idx = 0; idx < lhs.size(); idx++) {
            result.push_back(lhs[idx] * rhs[idx]);
        }
        return result;
    }

    vector<double> inverse( const vector<double> &initial);

    double inverse( double initial);
}

#endif //MAIN_CPP_VECTOR_HPP
