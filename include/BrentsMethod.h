#pragma once

#include <functional>

namespace BrentsMethod {

/// Function to bracket a root of a given function. Adapted from Numerical Recipes in C
///
/// @param f Reference to function to find bracketing interval
/// @param[out] x1 Reference one bound
/// @param[out] x2 Reference to other bound
void bracket(std::function<double(double)> const & f, double & x1, double & x2);

/// Finds the root of a function using Brent's method. Adapted from Numerical Recipes in C
///
/// @param f Reference to function to find root of
/// @param x1 One end of bracketing interval
/// @param x2 Other end of bracketing interval
/// @param tolerance Root finding tolerance (default is 1e-12)
double root(std::function<double(double)> const & f, double x1, double x2, double tol = 1.0e-12);

} // namespace BrentsMethod
