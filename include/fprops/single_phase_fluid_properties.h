// SPDX-FileCopyrightText: 2022 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#pragma once

#include "fprops/state.h"

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
    [[nodiscard]] State rho_T(double rho, double T) const;

    /// Compute thermodynamical state given density and pressure
    ///
    /// @param rho Density \f$[kg/m^3]\f$
    /// @param p Pressure \f$[Pa]\f$
    [[nodiscard]] State rho_p(double rho, double p) const;

    /// Compute thermodynamical state given pressure and temperature
    ///
    /// @param p Pressure \f$[Pa]\f$
    /// @param T Temperature \f$[K]\f$
    [[nodiscard]] State p_T(double p, double T) const;

    /// Compute thermodynamical state given specific volume and internal energy
    ///
    /// @param v Specific volume \f$[m^3/kg]\f$
    /// @param u Specific internal energy \f$[J/kg]\f$
    [[nodiscard]] State v_u(double v, double u) const;

    /// Compute thermodynamical state given specific enthalpy and entropy
    ///
    /// @param h Specific enthalpy \f$[J/kg]\f$
    /// @param s Entropy \f$[J/(kg-K)]\f$
    [[nodiscard]] State h_s(double h, double s) const;
};

} // namespace fprops
