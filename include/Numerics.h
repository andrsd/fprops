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

namespace newton {

/// Finds the root of a function using Newton's method
///
/// @param x0 Initial guess
/// @param f Function to find root of
/// @param df Derivative of `f`
/// @param tolerance Root finding tolerance (default is 1e-12)
double root(double x0,
            std::function<double(double)> const & f,
            std::function<double(double)> const & df,
            double tol = 1.0e-12);

} // namespace newton

} // namespace fprops
