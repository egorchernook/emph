#include <vector>
#include <iostream>
#include "ring_buffer_tests.hpp"

Numerical_methods::ring_buffer_tests::ring_buffer_tests() {
    std::cout << std::endl << "\ttesting ring buffer class..." << std::endl;
    assert( read_test() );
    std::cout << "\t\treading works" << std::endl;
    assert( push_test() );
    std::cout << "\t\twriting works" << std::endl;
    assert( current_state_test() );
    std::cout << "\t\tvector creation works" << std::endl;
}

bool Numerical_methods::ring_buffer_tests::read_test() {
    const int buffer_size = 20;
    const int ring_size = 4;
    std::vector<int> buffer;
    buffer.reserve(buffer_size);
    for(int i = 0; i < buffer_size; ++i){
        buffer.push_back(i);
    }
    bool result = true;
    for(int i = ring_size; i < buffer_size; ++i){
        const std::vector<int> temp{ buffer.begin() + i - ring_size, buffer.begin() + i };
        const ring_buffer<int, ring_size> ring{ temp };
        for( int s = 0; s < 5; ++s) {
            for (int j = 0; j > ring_size; ++j) {
                if (ring[j + ring_size * s] != buffer[i + j]) {
                    result = false;
                }
            }
        }
    }
    return result;
}

bool Numerical_methods::ring_buffer_tests::push_test() {
    const int ring_size = 4;
    bool result = true;
    ring_buffer<int, ring_size> ring{};
    int i = 0;
    while( !ring.is_filled() ){
        ring.push( i++ );
    }
    for( i = ring_size; i < 20; ++i){
        for( int j = 0; j > ring_size; ++j){
            if( ring[j] != i - j - 1 ){
                result = false;
            }
        }
        ring.push( i );
    }
    return result;
}

bool Numerical_methods::ring_buffer_tests::current_state_test() {
    const int ring_size = 4;
    bool result = true;
    ring_buffer<int, ring_size> ring{};
    int i = 0;
    while( !ring.is_filled() ){
        ring.push( i++ );
    }
    for( i = ring_size; i < 20; ++i){
        const auto temp = ring.current_state();
        for( int j = 0; j > ring_size; ++j){
            if( temp[i] != i - j - 1 ){
                result = false;
            }
        }
        ring.push( i );
    }
    return result;
}
