#include "Utils.h"
#include <cassert>
#include <cmath>

namespace fprops {
namespace utils {

double
normalize_interval_location(double lo, double hi, double x)
{
    double d = hi - lo;
    double xs = x - lo;
    if (std::abs(d) > 1e-10)
        xs /= d;
    return xs;
}

} // namespace utils
} // namespace fprops
