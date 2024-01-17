#pragma once

#include <functional>

namespace fprops {

inline double
sqr(double x)
{
    return x * x;
}

inline double
cb(double x)
{
    return x * x * x;
}

namespace math {

template <typename T>
T
pow(T x, int e)
{
    bool neg = false;
    T result = 1.0;

    if (e < 0) {
        neg = true;
        e = -e;
    }

    while (e) {
        // if bit 0 is set multiply the current power of two factor of the exponent
        if (e & 1)
            result *= x;
        // x is incrementally set to consecutive powers of powers of two
        x *= x;
        // bit shift the exponent down
        e >>= 1;
    }

    return neg ? 1.0 / result : result;
}

} // namespace math

namespace newton {

/// Finds the root of a function using Newton's method
///
/// @param x0 Initial guess
/// @param f Function to find root of
/// @param df Derivative of `f`
/// @param tol Root finding tolerance (default is 1e-12)
double root(double x0,
            std::function<double(double)> const & f,
            std::function<double(double)> const & df,
            double tol = 1.0e-12);

} // namespace newton

} // namespace fprops
