#ifndef MAIN_CPP_RING_BUFFER_HPP
#define MAIN_CPP_RING_BUFFER_HPP

#include <cstdlib>
#include "vector.hpp"
#include "solution_t_concept.hpp"
namespace Numerical_methods {

    template<
            typename solution_t,
            std::size_t elements_amount
    >
    class ring_buffer {
        using  base_container_t = std::vector< solution_t>;
    private:
        base_container_t data;
        std::size_t start = 0;
        bool full = false;
    public:

        ring_buffer() : data() {}

        ring_buffer(const base_container_t &_container, std::size_t _start = 0, bool _full = false)
                : data(_container),
                  start(_start),
                  full(_full) {
            assert(_container.size() == elements_amount);
        }

        solution_t operator[](std::size_t idx) const noexcept {
            assert( is_filled() );
            return data[(start + idx) % elements_amount];
        }

        solution_t at(std::size_t idx) const noexcept {
            assert( is_filled() );
            return data[(start + idx) % elements_amount];
        }

        void push(const solution_t &value) {
            if(!full){
                data.push_back(value);
            } else {
                data.at(start) = value;
            }
            if ((start + 1) / elements_amount != 0) {
                full = true;
            }
            ++start %= elements_amount;
        }

        std::vector<solution_t> current_state(){
            //assert( is_filled());
            std::vector<solution_t> result;
            if(is_filled() ) {
                result.reserve(elements_amount);
                for (int i = 0; i < elements_amount; ++i) {
                    result.push_back(this->at(i));
                }
            } else {
                for( int i = 0; i < elements_amount; ++i){
                    result.push_back(data[i]);
                }
            }
            return result;
        }

        [[nodiscard]] solution_t first() noexcept {
            return data[(start - 1) % elements_amount];
        }

        [[nodiscard]] solution_t last() noexcept {
            if( full ) {
                return data[start];
            } else {
                return data[0];
            }
        }

        [[nodiscard]] bool is_filled() const noexcept {
            return full;
        }

        void clear(){
            data.clear();
            start = 0;
            full = false;
        }
    };
}

#endif //MAIN_CPP_RING_BUFFER_HPP
