// SPDX-FileCopyrightText: 2024 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#pragma once

#include "fprops/single_phase_fluid_properties.h"
#include <string>

namespace fprops {
/// Construct fluid properties from a fluid name
///
/// @param name Fluid name
/// @return Fluid properties
SinglePhaseFluidProperties * from_name(const std::string & name);

} // namespace fprops
