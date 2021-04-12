#ifndef MAIN_CPP_MIN_AND_MAX_FUNCTIONS_HPP
#define MAIN_CPP_MIN_AND_MAX_FUNCTIONS_HPP
#include <cstdlib>
#include <type_traits>

namespace Numerical_methods {

    double max(double arg);
    double min(double arg);
    double abs(double arg);

    template<typename T>
    concept sequence_container = requires( T arg){
        arg[0];
        arg.at(0);
        arg.size();
    };
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
    Class abs(const Class &arg) {
        Class result = arg;
        for (int i = 0; i < arg.size(); ++i) {
            result[i] = abs(arg[i]);
        }
        return result;
    }
}

#endif //MAIN_CPP_MIN_AND_MAX_FUNCTIONS_HPP
