// SPDX-FileCopyrightText: 2022 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#include "fprops/single_phase_fluid_properties.h"
#include "fprops/exception.h"

namespace fprops {

SinglePhaseFluidProperties::SinglePhaseFluidProperties() {}

State
SinglePhaseFluidProperties::rho_T(double rho, double T) const
{
    throw Exception("Calling generic fluid properties");
}

State
SinglePhaseFluidProperties::rho_p(double rho, double p) const
{
    throw Exception("Calling generic fluid properties");
}

State
SinglePhaseFluidProperties::p_T(double p, double T) const
{
    throw Exception("Calling generic fluid properties");
}

State
SinglePhaseFluidProperties::v_u(double v, double u) const
{
    throw Exception("Calling generic fluid properties");
}

State
SinglePhaseFluidProperties::h_s(double h, double s) const
{
    throw Exception("Calling generic fluid properties");
}

} // namespace fprops
