#pragma once

#include "SinglePhaseFluidProperties.h"

namespace fprops {

/// Base class for fluid properties based on Helmholtz equation of state
///
/// This class is based on `HelmholtzFluidProperties.h` from `idaholab/moose/fluid_properties`
/// module
class Helmholtz : public SinglePhaseFluidProperties {
public:
    /// @param R Universal gas constant \f$[J/(mol-K)]\f$
    /// @param M Molar mass \f$[kg/mol]\f$
    /// @param rho_c Critical density \f$[kg/m^3]\f$
    /// @param T_c Critical temperature \f$[K]\f$
    Helmholtz(double R, double M, double rho_c, double T_c);

    [[nodiscard]] Props p_T(double p, double T) const override;
    [[nodiscard]] Props v_u(double v, double u) const override;

protected:
    /// Helmholtz free energy
    ///
    /// @param delta Scaled density \f$[-]\f$
    /// @param tau Scaled temperature \f$[-]\f$
    /// @return Helmholtz free energy (\f$\alpha\f$)
    [[nodiscard]] virtual double alpha(double delta, double tau) const = 0;

    /// Derivative of Helmholtz free energy wrt delta
    ///
    /// @param delta Scaled density \f$[-]\f$
    /// @param tau Scaled temperature \f$[-]\f$
    /// @return Derivative of Helmholtz free energy wrt delta
    [[nodiscard]] virtual double dalpha_ddelta(double delta, double tau) const = 0;

    /// Derivative of Helmholtz free energy wrt tau
    ///
    /// @param delta Scaled density \f$[-]\f$
    /// @param tau Scaled temperature \f$[-]\f$
    /// @return Derivative of Helmholtz free energy wrt tau
    [[nodiscard]] virtual double dalpha_dtau(double delta, double tau) const = 0;

    /// Second derivative of Helmholtz free energy wrt delta
    ///
    /// @param delta Scaled density \f$[-]\f$
    /// @param tau Scaled temperature \f$[-]\f$
    /// @return Second derivative of Helmholtz free energy wrt delta
    [[nodiscard]] virtual double d2alpha_ddelta2(double delta, double tau) const = 0;

    /// Second derivative of Helmholtz free energy wrt tau
    ///
    /// @param delta Scaled density \f$[-]\f$
    /// @param tau Scaled temperature \f$[-]\f$
    /// @return Second derivative of Helmholtz free energy wrt tau
    [[nodiscard]] virtual double d2alpha_dtau2(double delta, double tau) const = 0;

    /// Second derivative of Helmholtz free energy wrt delta and tau
    ///
    /// @param delta Scaled density \f$[-]\f$
    /// @param tau Scaled temperature \f$[-]\f$
    /// @return Second derivative of Helmholtz free energy wrt delta and tau
    [[nodiscard]] virtual double d2alpha_ddeltatau(double delta, double tau) const = 0;

    /// Density given pressure and temperature
    ///
    /// @param p Pressure \f$[Pa]\f$
    /// @param T Temperature \f$[K]\f$
    /// @return Density \f$[kg/m^3]\f$
    [[nodiscard]] double rho_from_p_T(double p, double T) const;

    /// Tau given specific volume and internal energy
    ///
    /// @param v Specific volume \f$[m^3/kg]\f$
    /// @param u Specific internal energy \f$[J/kg]\f$
    /// @return Tau (scaled temperature)
    [[nodiscard]] double tau_from_v_u(double v, double u) const;

    /// Dynamic viscosity
    ///
    /// @param rho Density \f$[kg/m^3]\f$
    /// @param T Temperature \f$[K]\f$
    /// @return Dynamic viscosity \f$[Pa-s]\f$
    [[nodiscard]] virtual double mu_from_rho_T(double rho, double T) const = 0;

    /// Thermal conductivity
    ///
    /// @param rho Density \f$[kg/m^3]\f$
    /// @param T Temperature \f$[K]\f$
    /// @return Thermal conductivity \f$[W/(m-K)]\f$
    [[nodiscard]] virtual double k_from_rho_T(double rho, double T) const = 0;

    /// Universal gas constant \f$[J/(mol-K)]\f$
    const double R;
    /// Molar mass \f$[kg/mol]\f$
    const double M;
    /// Critical density \f$[kg/m^3]\f$
    const double rho_c;
    /// Critical temperature \f$[K]\f$
    const double T_c;
};

} // namespace fprops
