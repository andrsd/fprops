// SPDX-FileCopyrightText: 2022 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#pragma once

#include "fprops/single_phase_fluid_properties.h"

namespace fprops {

/// Ideal gas
///
class IdealGas {
public:
    /// Constructor
    ///
    /// @param gamma Adiabatic index (ratio of specific heats cp/cv)
    /// @param molar_mass Molar mass \f$[kg/mol]\f$
    IdealGas(double gamma, double molar_mass);

    /// Get adiabatic index
    ///
    /// @return Adiabatic index (ratio of specific heats cp/cv)
    double gamma() const;

    /// Get molar mass
    ///
    /// @return Molar mass \f$[kg/mol]\f$
    double molar_mass() const;

    /// Get specific gas constant
    ///
    /// @return Specific gas constant
    double R_specific() const;

    /// Get specific heat at constant pressure
    ///
    /// @return Specific heat at constant pressure \f$[J/(kg-K)]\f$
    double cp() const;

    /// Specific heat at constant volume \f$[J/(kg-K)]\f$
    ///
    /// @return Specific heat at constant volume \f$[J/(kg-K)]\f$
    double cv() const;

    /// Get dynamic viscosity
    ///
    /// @return Dynamic viscosity \f$[Pa-s]\f$
    double mu() const;

    /// Get thermal conductivity
    ///
    /// @return Thermal conductivity \f$[W/(m-K)]\f$
    double k() const;

    /// Set dynamic viscosity
    ///
    /// @param mu Dynamic viscosity \f$[Pa-s]\f$
    void set_mu(double mu);

    /// Set thermal conductivity
    ///
    /// @param k Thermal conductivity \f$[W/(m-K)]\f$
    void set_k(double k);

    [[nodiscard]] State rho_T(double rho, double T) const;
    [[nodiscard]] State rho_p(double rho, double p) const;
    [[nodiscard]] State p_T(double p, double T) const;
    [[nodiscard]] State v_u(double v, double u) const;
    [[nodiscard]] State h_s(double h, double s) const;
    [[nodiscard]] State v_h(double v, double h) const;

private:
    /// Adiabatic index (ratio of specific heats cp/cv)
    double m_gamma;
    /// Molar mass \f$[kg/mol]\f$
    double m_molar_mass;
    /// Specific gas constant (R / molar mass)
    double m_R_specific;
    /// Specific heat at constant pressure \f$[J/(kg-K)]\f$
    double m_cp;
    /// Specific heat at constant volume \f$[J/(kg-K)]\f$
    double m_cv;
    /// Dynamic viscosity \f$[Pa-s]\f$
    double m_mu;
    /// Thermal conductivity \f$[W/(m-K)]\f$
    double m_k;
};

} // namespace fprops
