#include <iostream>
#include "matrix_tests.hpp"

Numerical_methods::matrix_tests::matrix_tests() {
    std::cout << std::endl << "\ttesting matrix class..." << std::endl;
    assert( creation_tests() );
    std::cout << "\t\tconstruction works" << std::endl;
    assert( two_parameters_at_test() );
    std::cout << "\t\tat with two parameters works" << std::endl;
    assert( direct_assignment_test() );
    std::cout << "\t\tdirect assignment works" << std::endl;
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
    std::cout << "\t\tdivision by a number assignment works" << std::endl;
    assert( two_matrix_multiplication_test() );
    std::cout << "\t\tmultiplication of two matrices works" << std::endl;
    assert( inverse_Gauss_test() );
    std::cout << "\t\tGauss inversion works" << std::endl;
}

bool Numerical_methods::matrix_tests::creation_tests() {
    constexpr int size = 3;

    auto** matrix = new(std::nothrow) double * [size] ;
    for( int i = 0; i < size; ++i){
        matrix[i] = new(std::nothrow) double [size];
    }
    matrix[0][0] = 2.58; matrix[0][1] = 2.93; matrix[0][2] = 3.13;
    matrix[1][0] = 1.32; matrix[1][1] = 1.55; matrix[1][2] = 1.58;
    matrix[2][0] = 2.09; matrix[2][1] = 2.25; matrix[2][2] = 2.34;

    const Numerical_methods::matrix<double, size, size> a = { { 2.58 , 2.93 , 3.13 },
                                                              { 1.32 , 1.55 , 1.58 },
                                                              { 2.09 , 2.25 , 2.34 }  };
    Numerical_methods::matrix<double, size, size> b;
    b[0][0] = 2.58; b[0][1] = 2.93; b[0][2] = 3.13;
    b[1][0] = 1.32; b[1][1] = 1.55; b[1][2] = 1.58;
    b[2][0] = 2.09; b[2][1] = 2.25; b[2][2] = 2.34;

    auto c = a;
    bool result = true;
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            const bool values_differ =
                    std::abs(matrix[i][j] - a[i][j]) > std::numeric_limits<double>::epsilon() &&
                    std::abs(a[i][j] - b[i][j])  > std::numeric_limits<double>::epsilon() &&
                    std::abs(b[i][j] - c[i][j])  > std::numeric_limits<double>::epsilon() &&
                    std::abs(a[i][j] - c[i][j])  > std::numeric_limits<double>::epsilon();
            if( values_differ ){
                result = false;
            }
        }
    }
    for (int i = 0; i < size; i++)
        delete[] matrix[i];
    delete[] matrix;
    return result;
}

bool Numerical_methods::matrix_tests::two_parameters_at_test() {
    constexpr int size = 3;
    const Numerical_methods::matrix<double, size, size> a{ 5.0 };
    bool result = true;
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            if( std::abs(a.at(i,j) - 5.0) > std::numeric_limits<double>::epsilon() * 1'000) {
                result = false;
            }
        }
    }
    return result;
}

bool Numerical_methods::matrix_tests::direct_assignment_test() {
    constexpr int size = 3;
    const Numerical_methods::matrix<double, size, size> a{ 5.0 };
    Numerical_methods::matrix<double, size, size> b{0.0};
    b = a;
    bool result = true;
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            if( std::abs(b[i][j] - 5.0) > std::numeric_limits<double>::epsilon() * 1'000) {
                result = false;
            }
        }
    }
    return result;
}

bool Numerical_methods::matrix_tests::addition_assignment_test() {
    constexpr int size = 3;
    Numerical_methods::matrix<double, size, size> a{ 5.0 };
    a += {{1.0,1.0,1.0},{1.0,1.0,1.0},{1.0,1.0,1.0}};
    bool result = true;
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            if( std::abs(a[i][j] - 6.0) > std::numeric_limits<double>::epsilon() * 1'000) {
                result = false;
            }
        }
    }
    return result;
}

bool Numerical_methods::matrix_tests::subtraction_assignment_test() {
    constexpr int size = 3;
    Numerical_methods::matrix<double, size, size> a{ 5.0 };
    a -= {{1.0,1.0,1.0},{1.0,1.0,1.0},{1.0,1.0,1.0}};
    bool result = true;
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            if( std::abs(a[i][j] - 4.0) > std::numeric_limits<double>::epsilon() * 1'000) {
                result = false;
            }
        }
    }
    return result;
}

bool Numerical_methods::matrix_tests::multiplication_with_number_assignment_test() {
    constexpr int size = 3;
    Numerical_methods::matrix<double, size, size> a{ 5.0 };
    a *= 2.0;
    bool result = true;
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            if( std::abs(a[i][j] - 10.0) > std::numeric_limits<double>::epsilon() * 1'000) {
                result = false;
            }
        }
    }
    return result;
}

bool Numerical_methods::matrix_tests::division_by_number_assignment_test() {
    constexpr int size = 3;
    Numerical_methods::matrix<double, size, size> a{ 10.0 };
    a /= 5.0;
    bool result = true;
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            if( std::abs(a[i][j] - 2.0) > std::numeric_limits<double>::epsilon() * 1'000) {
                result = false;
            }
        }
    }
    return result;
}

bool Numerical_methods::matrix_tests::two_matrix_multiplication_test() {
    constexpr int height = 2;
    constexpr int width = 3;
    const matrix<double, height, width> lhs = { { 1.0, 2.0, 3.0 },
                                                { 1.0, 2.0, 3.0} };
    const matrix<double, width, height> rhs = { { 1.0, 2.0} ,
                                                { 1.0, 2.0} ,
                                                { 1.0, 2.0} };
    const matrix<double, height, height> right_answer = { { 6.0, 12.0 },
                                                          { 6.0, 12.0 } };
    const auto answer = lhs * rhs;
    bool result = true;
    for(int i = 0; i < height; ++i){
        for(int j = 0; j < height; ++j){
            if( std::abs( right_answer[i][j] - answer[i][j] ) > std::numeric_limits<double>::epsilon() * 1'000 )
            {
                result = false;
            }
        }
    }
    return result;
}

bool Numerical_methods::matrix_tests::inverse_Gauss_test() {
    constexpr int size = 3;
    constexpr double error = 0.001;
    const Numerical_methods::matrix<double, size, size> a = { { 2.58, 2.93, 3.13 } ,
                                                              { 1.32, 1.55, 1.58 } ,
                                                              { 2.09, 2.25, 2.34 } };
    auto c = inverse_Gauss(a);

    auto B = a * c;
    bool result = true;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if ( i == j ) {
                if ( std::abs(B[i][j] - 1.0) > error ){
                    result = false;
                }
            } else {
                if ( std::abs( B[i][j] ) > error ) {
                    result = false;
                }
            }
        }
    }
    return result;
}

bool Numerical_methods::matrix_tests::addition_test() {
    constexpr int size = 3;
    Numerical_methods::matrix<double, size, size> a{ 1.0 };
    Numerical_methods::matrix<double, size, size> b{ 1.0 };
    auto c = a + b;
    bool result = true;
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            if( std::abs( c[i][j] - 2.0 ) > std::numeric_limits<double>::epsilon() * 1'000 )
            {
                result = false;
            }
        }
    }
    return result;
}

bool Numerical_methods::matrix_tests::subtraction_test() {
    constexpr int size = 3;
    Numerical_methods::matrix<double, size, size> a{ 3.0 };
    Numerical_methods::matrix<double, size, size> b{ 1.0 };
    auto c = a - b;
    bool result = true;
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            if( std::abs( c[i][j] - 2.0 ) > std::numeric_limits<double>::epsilon() * 1'000 )
            {
                result = false;
            }
        }
    }
    return result;
}

bool Numerical_methods::matrix_tests::multiplication_with_number_test() {
    constexpr int size = 3;
    Numerical_methods::matrix<double, size, size> a{ 1.0 };
    auto b = a * 2.0;
    auto c = 3.0 * a;
    bool result = true;
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            if( std::abs( c[i][j] - 3.0 ) > std::numeric_limits<double>::epsilon() * 1'000 &&
                std::abs( b[i][j] - 2.0 ) > std::numeric_limits<double>::epsilon() * 1'000 )
            {
                result = false;
            }
        }
    }
    return result;
}

bool Numerical_methods::matrix_tests::division_by_number_test() {
    constexpr int size = 3;
    Numerical_methods::matrix<double, size, size> a{ 6.0 };
    auto b = a / 2.0;
    bool result = true;
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            if( std::abs( b[i][j] - 3.0 ) > std::numeric_limits<double>::epsilon() * 1'000 )
            {
                result = false;
            }
        }
    }
    return result;
}

