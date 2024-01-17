#pragma once

#include <vector>

namespace fprops {
namespace utils {

/// Find interval index for given location
///
/// @param range Array (in ascending order) of locations
/// @param x Location we re looking for
/// @return Index of the interval in the `range` array
///
/// Example:
/// - range = [1, 2.5, 5]
/// - x = 1.2
/// - return value 0 (1.2 is on the 0th interval)
template <typename T>
std::size_t
interval_index(const T & range, double x)
{
    assert(range.size() >= 2);

    std::size_t lo = 0;
    std::size_t hi = range.size() - 1;
    while (hi - lo > 1) {
        unsigned int mid = (hi + lo) / 2;
        if (x < range[mid])
            hi = mid;
        else
            lo = mid;
    }
    return lo;
}

///
double normalize_interval_location(double lo, double hi, double x);

} // namespace utils
} // namespace fprops
