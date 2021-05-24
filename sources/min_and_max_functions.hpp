#ifndef MAIN_CPP_MIN_AND_MAX_FUNCTIONS_HPP
#define MAIN_CPP_MIN_AND_MAX_FUNCTIONS_HPP
#include <cstdlib>
#include <type_traits>

namespace Numerical_methods {

    double max(double arg);
    double min(double arg);
    double lowest(double arg);
    double abs(double arg);
    std::size_t size( double);

    #if defined __GNUC__ && __GNUC__ < 7
        #define sequence_container typename
    #else
        template<typename T>
        concept sequence_container = requires( T arg){
            arg[0];
            arg.at(0);
            arg.size();
        };
    #endif

    template<sequence_container Class>
    double max(const Class &arg) {
        double result = max(arg[0]);
        for (int i = 0; i < arg.size(); ++i) {
            if (result < max(arg[i])) {
                result = max(arg[i]);
            }
        }
        return result;
    }

    template<sequence_container Class>
    double min(const Class &arg) {
        double result = min(arg[0]);
        for (int i = 0; i < arg.size(); ++i) {
            if (result > min(arg[i])) {
                result = min(arg[i]);
            }
        }
        return result;
    }

    template<sequence_container Class>
    double lowest(const Class &arg) {
        double result = lowest(arg[0]);
        for (int i = 0; i < arg.size(); ++i){
            if( abs(result) > abs(lowest(arg[i]))){
                result = lowest(arg[i]);
            }
        }
        return result;
    }

    template<sequence_container Class>
    Class abs(const Class &arg) {
        Class result = arg;
        for (int i = 0; i < arg.size(); ++i) {
            result[i] = abs(arg[i]);
        }
        return result;
    }

    template<sequence_container Class>
    std::size_t size( const Class &arg){
        return arg.size();
    }
}

#endif //MAIN_CPP_MIN_AND_MAX_FUNCTIONS_HPP
