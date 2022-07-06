#pragma once

#include "SinglePhaseFluidProperties.h"

namespace fprops {

/// Ideal gas
///
class IdealGas : public SinglePhaseFluidProperties {
public:
    IdealGas(double gamma, double molar_mass);

    /// Set dynamic viscosity
    ///
    /// @param mu Dynamic viscosity
    void set_mu(double mu);

    /// Set thermal conductivity
    ///
    /// @param k Thermal conductivity
    void set_k(double k);

    virtual Props p_T(double p, double T) override;
    virtual Props v_u(double v, double u) override;

protected:
    /// Adiabatic index (ratio of specific heats cp/cv)
    double gamma;
    /// Molar mass
    double molar_mass;
    /// Specific gas constant (R / molar mass)
    double R_specific;
    /// Specific heat at constant pressure
    double cp;
    /// Specific heat at constant volume
    double cv;
    /// Dynamic viscosity
    double mu;
    /// Thermal conductivity
    double k;
};

} // namespace fprops
