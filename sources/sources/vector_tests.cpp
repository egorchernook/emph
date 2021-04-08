#include <iostream>
#include <cassert>
#include "vector_tests.hpp"

Numerical_methods::vector_tests::vector_tests() {
    std::cout << std::endl << "\ttesting vector class..." << std::endl;
    assert( addition_assignment_test() );
    std::cout << "\t\taddition assignment works" << std::endl;
    assert( subtraction_assignment_test() );
    std::cout << "\t\tsubtraction assignment works" << std::endl;
    assert( multiplication_with_number_assignment_test() );
    std::cout << "\t\tmultiplication on a number assignment works" << std::endl;
    assert( division_by_number_assignment_test() );
    std::cout << "\t\tdivision by a number assignment works" << std::endl;
    assert( addition_test() );
    std::cout << "\t\taddition works" << std::endl;
    assert( subtraction_test() );
    std::cout << "\t\tsubtraction works" << std::endl;
    assert( multiplication_with_number_test() );
    std::cout << "\t\tmultiplication with a number works" << std::endl;
    assert( division_by_number_test() );
}

bool Numerical_methods::vector_tests::addition_assignment_test() {
    const Numerical_methods::vector<double> a = { 1.0, 2.0, 3.0 };
    Numerical_methods::vector<double> b = { 3.0, 2.0, 1.0 };
    b += a;
    bool result = true;
    for(int i = 0; i < a.size(); ++i){
        if( std::abs( b[i] - 4.0) > std::numeric_limits<double>::epsilon() * 1'000) {
            result = false;
        }
    }
    return result;
}

bool Numerical_methods::vector_tests::subtraction_assignment_test() {
    const Numerical_methods::vector<double> a = { 1.0, 1.0, 1.0 };
    Numerical_methods::vector<double> b = { 3.0, 3.0, 3.0 };
    b -= a;
    bool result = true;
    for(int i = 0; i < a.size(); ++i){
        if( std::abs( b[i] - 2.0) > std::numeric_limits<double>::epsilon() * 1'000) {
            result = false;
        }
    }
    return result;
}

bool Numerical_methods::vector_tests::multiplication_with_number_assignment_test() {
    Numerical_methods::vector<double> b = { 3.0, 3.0, 3.0 };
    b *= 2.0;
    bool result = true;
    for(int i = 0; i < b.size(); ++i){
        if( std::abs( b[i] - 6.0) > std::numeric_limits<double>::epsilon() * 1'000) {
            result = false;
        }
    }
    return result;
}

bool Numerical_methods::vector_tests::division_by_number_assignment_test() {
    Numerical_methods::vector<double> b = { 4.0, 4.0, 4.0 };
    b /= 2.0;
    bool result = true;
    for(int i = 0; i < b.size(); ++i){
        if( std::abs( b[i] - 2.0) > std::numeric_limits<double>::epsilon() * 1'000) {
            result = false;
        }
    }
    return result;
}

bool Numerical_methods::vector_tests::addition_test() {
    const Numerical_methods::vector<double> a = { 1.0, 2.0, 3.0 };
    const Numerical_methods::vector<double> b = { 3.0, 2.0, 1.0 };
    auto c = a + b;
    bool result = true;
    for(int i = 0; i < a.size(); ++i){
        if( std::abs( c[i] - 4.0) > std::numeric_limits<double>::epsilon() * 1'000) {
            result = false;
        }
    }
    return result;
}

bool Numerical_methods::vector_tests::subtraction_test() {
    const Numerical_methods::vector<double> a = { 3.0, 3.0, 3.0 };
    const Numerical_methods::vector<double> b = { 1.0, 1.0, 1.0 };
    auto c = a - b;
    bool result = true;
    for(int i = 0; i < a.size(); ++i){
        if( std::abs( c[i] - 2.0) > std::numeric_limits<double>::epsilon() * 1'000) {
            result = false;
        }
    }
    return result;
}

bool Numerical_methods::vector_tests::multiplication_with_number_test() {
    const Numerical_methods::vector<double> a = { 3.0, 3.0, 3.0 };
    auto b = 3.0 * a;
    auto c = a * 2.0;
    bool result = true;
    for(int i = 0; i < a.size(); ++i){
        if( std::abs( b[i] - 9.0) > std::numeric_limits<double>::epsilon() * 1'000 &&
            std::abs( c[i] - 6.0) > std::numeric_limits<double>::epsilon() * 1'000 ) {
            result = false;
        }
    }
    return result;
}

bool Numerical_methods::vector_tests::division_by_number_test() {
    const Numerical_methods::vector<double> a = { 4.0, 4.0, 4.0 };
    auto c = a / 2.0;
    bool result = true;
    for(int i = 0; i < a.size(); ++i){
        if( std::abs( c[i] - 2.0) > std::numeric_limits<double>::epsilon() * 1'000 ) {
            result = false;
        }
    }
    return result;
}
