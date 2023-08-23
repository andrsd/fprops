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

    /// Set dynamic viscosity
    ///
    /// @param mu Dynamic viscosity \f$[Pa-s]\f$
    void set_mu(double mu);

    /// Set thermal conductivity
    ///
    /// @param k Thermal conductivity \f$[W/(m-K)]\f$
    void set_k(double k);

    [[nodiscard]] Props rho_T(double rho, double T) const override;
    [[nodiscard]] Props rho_p(double rho, double p) const override;
    [[nodiscard]] Props p_T(double p, double T) const override;
    [[nodiscard]] Props v_u(double v, double u) const override;
    [[nodiscard]] Props h_s(double h, double s) const override;

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
