// SPDX-FileCopyrightText: 2022 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#pragma once

#include "fprops/State.h"

namespace fprops {

/// Base class for single phase fluids
///
class SinglePhaseFluidProperties {
public:
    SinglePhaseFluidProperties();
    virtual ~SinglePhaseFluidProperties() = default;

    /// Compute thermodynamical state given density and temperature
    ///
    /// @param rho Density \f$[kg/m^3]\f$
    /// @param T Temperature \f$[K]\f$
    [[nodiscard]] virtual State rho_T(double rho, double T) const = 0;

    /// Compute thermodynamical state given density and pressure
    ///
    /// @param rho Density \f$[kg/m^3]\f$
    /// @param p Pressure \f$[Pa]\f$
    [[nodiscard]] virtual State rho_p(double rho, double p) const = 0;

    /// Compute thermodynamical state given pressure and temperature
    ///
    /// @param p Pressure \f$[Pa]\f$
    /// @param T Temperature \f$[K]\f$
    [[nodiscard]] virtual State p_T(double p, double T) const = 0;

    /// Compute thermodynamical state given specific volume and internal energy
    ///
    /// @param v Specific volume \f$[m^3/kg]\f$
    /// @param u Specific internal energy \f$[J/kg]\f$
    [[nodiscard]] virtual State v_u(double v, double u) const = 0;

    /// Compute thermodynamical state given specific enthalpy and entropy
    ///
    /// @param h Specific enthalpy \f$[J/kg]\f$
    /// @param s Entropy \f$[J/(kg-K)]\f$
    [[nodiscard]] virtual State h_s(double h, double s) const = 0;
};

} // namespace fprops
