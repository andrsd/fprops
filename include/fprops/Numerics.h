// SPDX-FileCopyrightText: 2022 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#pragma once

#include <functional>
#include <stdexcept>
#include <type_traits>

namespace fprops {

namespace math {

template <int N, typename T>
struct pow_impl {
    static inline T
    value(T x)
    {
        if (N % 2)
            return x * pow_impl<N - 1, T>::value(x);
        T x_n_half = pow_impl<N / 2, T>::value(x);
        return x_n_half * x_n_half;
    }
};

template <typename T>
struct pow_impl<1, T> {
    static inline T
    value(T x)
    {
        return x;
    }
};

template <typename T>
struct pow_impl<0, T> {
    static inline T
    value(T)
    {
        return 1;
    }
};

template <int N, typename T>
inline T
pow(T x)
{
    return pow_impl<N, T>::value(x);
}

/// Compute `e` power of `x`
///
/// @param x Operand
/// @param e Exponent (integer)
/// @return `x` to the power of `e`
template <typename T>
inline T
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
