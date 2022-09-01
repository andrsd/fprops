#include "Numerics.h"
#include <cmath>
#include <stdexcept>

namespace fprops {

namespace newton {

/// The number of iterations during Newton's method until we declare non-convergence
static const unsigned int MAX_ITERATIONS = 50;
/// Numerical tolerance for the y_prime to prevent division by zero
static const double epsilon = 1e-10;

double
root(double x0,
     std::function<double(double)> const & f,
     std::function<double(double)> const & df,
     double tol)
{
    for (unsigned int iter = 0; iter < MAX_ITERATIONS; iter++) {
        double y = f(x0);
        double y_prime = df(x0);
        if (std::fabs(y_prime) < epsilon)
            break;

        double x1 = x0 - y / y_prime;
        if (std::fabs(x1 - x0) < tol)
            return x1;
        x0 = x1;
    }

    throw std::runtime_error("Newton's method failed to converge");
}

} // namespace newton

} // namespace fprops
