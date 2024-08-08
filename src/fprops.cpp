// SPDX-FileCopyrightText: 2024 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#include "fprops/fprops.h"
#include "fprops/air.h"
#include "fprops/carbon_dioxide.h"
#include "fprops/helium.h"
#include "fprops/nitrogen.h"
#include "fprops/exception.h"

namespace fprops {

std::string
to_lower(const std::string & name)
{
    std::string lower(name);
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    return lower;
}

SinglePhaseFluidProperties *
from_name(const std::string & name)
{
    auto n = to_lower(name);
    if (n == "air")
        return new Air();
    else if (n == "carbon_dioxide" || n == "co2")
        return new CarbonDioxide();
    else if (n == "helium" || n == "he")
        return new Helium();
    else if (n == "nitrogen" || n == "n2")
        return new Nitrogen();
    else
        throw Exception("Unknown fluid name '{}'", name);
}

} // namespace fprops
