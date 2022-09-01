#pragma once

#include "SinglePhaseFluidProperties.h"

namespace fprops {

/// Base class for fluid properties basedon Helmholtz equation of state
///
/// This class is based on HelmholtzFluidProperties.h from idaholab/moose/fluid_properties module
class Helmholtz : public SinglePhaseFluidProperties {
public:
    /// @param R Universal gas constant [J / (mol K)]
    /// @param M Molar mass [kg/mol]
    /// @param rho_c Critical density [kg/m^3]
    /// @param T_c Critical temperature [K]
    Helmholtz(double R, double M, double rho_c, double T_c);

    virtual Props p_T(double p, double T) override;
    virtual Props v_u(double v, double u) override;

protected:
    /// Helmholtz free energy
    ///
    /// @param delta Scaled density [-]
    /// @param tau Scaled temperature [-]
    /// @return alpha - Helmholtz free energy
    virtual double alpha(double delta, double tau) = 0;

    /// Derivative of Helmholtz free energy wrt delta
    ///
    /// @param delta Scaled density [-]
    /// @param tau Scaled temperature [-]
    /// @return Derivative of Helmholtz free energy wrt delta
    virtual double dalpha_ddelta(double delta, double tau) = 0;

    /// Derivative of Helmholtz free energy wrt tau
    ///
    /// @param delta Scaled density [-]
    /// @param tau Scaled temperature [-]
    /// @return Derivative of Helmholtz free energy wrt tau
    virtual double dalpha_dtau(double delta, double tau) = 0;

    /// Second derivative of Helmholtz free energy wrt delta
    ///
    /// @param delta Scaled density [-]
    /// @param tau Scaled temperature [-]
    /// @return Second derivative of Helmholtz free energy wrt delta
    virtual double d2alpha_ddelta2(double delta, double tau) = 0;

    /// Second derivative of Helmholtz free energy wrt tau
    ///
    /// @param delta Scaled density [-]
    /// @param tau Scaled temperature [-]
    /// @return Second derivative of Helmholtz free energy wrt tau
    virtual double d2alpha_dtau2(double delta, double tau) = 0;

    /// Second derivative of Helmholtz free energy wrt delta and tau
    ///
    /// @param delta Scaled density [-]
    /// @param tau Scaled temperature [-]
    /// @return Second derivative of Helmholtz free energy wrt delta and tau
    virtual double d2alpha_ddeltatau(double delta, double tau) = 0;

    double p_from_rho_T(double density, double temperature);

    double rho_from_p_T(double pressure, double temperature);

    /// Dynamic viscosity
    ///
    /// @param rho Density [kg/m^3]
    /// @param T Temperature [K]
    /// @return Dynamic viscosity [Pa.s]
    virtual double mu_from_rho_T(double rho, double T) = 0;

    /// Thermal conductivity
    ///
    /// @param rho Density [kg/m^3]
    /// @param T Temperature [K]
    /// @return Thermal conductivity [W/(m*K)]
    virtual double k_from_rho_T(double rho, double T) = 0;

    /// Universal gas constant [J / (mol K)]
    const double R;
    /// Molar mass [kg/mol]
    const double M;
    /// Critical density [kg/m^3]
    const double rho_c;
    /// Critical temperature [K]
    const double T_c;
};

} // namespace fprops
