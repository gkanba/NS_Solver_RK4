#ifndef __IG_H_UTIL_
#define __IG_H_UTIL_

namespace kanba
{
    template<typename T>
    struct Pair
    {
        T x;
        T y;
        Pair(const T x_i, const T y_i) : x(x_i), y(y_i) {}
    };
}

#endif
