// SPDX-FileCopyrightText: 2022 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#include "fprops/single_phase_fluid_properties.h"
#include "fprops/air.h"
#include "fprops/carbon_dioxide.h"
#include "fprops/helium.h"
#include "fprops/nitrogen.h"
#include "fprops/exception.h"

namespace fprops {

namespace {

std::string
to_lower(const std::string & name)
{
    std::string lower(name);
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    return lower;
}

} // namespace

SinglePhaseFluidProperties::~SinglePhaseFluidProperties()
{
    delete this->impl_;
}

State
SinglePhaseFluidProperties::rho_T(double rho, double T) const
{
    return this->impl_->rho_T(rho, T);
}

State
SinglePhaseFluidProperties::rho_p(double rho, double p) const
{
    return this->impl_->rho_p(rho, p);
}

State
SinglePhaseFluidProperties::p_T(double p, double T) const
{
    return this->impl_->p_T(p, T);
}

State
SinglePhaseFluidProperties::v_u(double v, double u) const
{
    return this->impl_->v_u(v, u);
}

State
SinglePhaseFluidProperties::h_s(double h, double s) const
{
    return this->impl_->h_s(h, s);
}

SinglePhaseFluidProperties
SinglePhaseFluidProperties::from_name(const std::string & name)
{
    auto n = to_lower(name);
    if (n == "air")
        return SinglePhaseFluidProperties(new Air());
    else if (n == "carbon_dioxide" || n == "co2")
        return SinglePhaseFluidProperties(new CarbonDioxide());
    else if (n == "helium" || n == "he")
        return SinglePhaseFluidProperties(new Helium());
    else if (n == "nitrogen" || n == "n2")
        return SinglePhaseFluidProperties(new Nitrogen());
    else
        throw Exception("Unknown fluid name '{}'", name);
}

} // namespace fprops
