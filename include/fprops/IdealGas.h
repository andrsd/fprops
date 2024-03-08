// SPDX-FileCopyrightText: 2022 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#pragma once

#include "SinglePhaseFluidProperties.h"

namespace fprops {

/// Ideal gas
///
class IdealGas : public SinglePhaseFluidProperties {
public:
    /// Constructor
    ///
    /// @param gamma Adiabatic index (ratio of specific heats cp/cv)
    /// @param molar_mass Molar mass \f$[kg/mol]\f$
    IdealGas(double gamma, double molar_mass);

    /// Get adiabatic index
    ///
    /// @return Adiabatic index (ratio of specific heats cp/cv)
    double get_gamma() const;

    /// Get specific gas constant
    ///
    /// @return Specific gas constant
    double get_specific_gas_constant() const;

    /// Set dynamic viscosity
    ///
    /// @param mu Dynamic viscosity \f$[Pa-s]\f$
    void set_mu(double mu);

    /// Set thermal conductivity
    ///
    /// @param k Thermal conductivity \f$[W/(m-K)]\f$
    void set_k(double k);

    [[nodiscard]] State rho_T(double rho, double T) const override;
    [[nodiscard]] State rho_p(double rho, double p) const override;
    [[nodiscard]] State p_T(double p, double T) const override;
    [[nodiscard]] State v_u(double v, double u) const override;
    [[nodiscard]] State h_s(double h, double s) const override;

protected:
    /// Adiabatic index (ratio of specific heats cp/cv)
    double gamma;
    /// Molar mass \f$[kg/mol]\f$
    double molar_mass;
    /// Specific gas constant (R / molar mass)
    double R_specific;
    /// Specific heat at constant pressure \f$[J/(kg-K)]\f$
    double cp;
    /// Specific heat at constant volume \f$[J/(kg-K)]\f$
    double cv;
    /// Dynamic viscosity \f$[Pa-s]\f$
    double mu;
    /// Thermal conductivity \f$[W/(m-K)]\f$
    double k;
};

} // namespace fprops
