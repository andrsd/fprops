// SPDX-FileCopyrightText: 2024 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#pragma once

#include "fprops/helmholtz.h"
#include "fprops/transport_models.h"

namespace fprops {

/// Oxygen fluid properties
///
/// References:
/// [1] R. Schmidt and W. Wagner. A New Form of the Equation of State for Pure Substances and its
///     Application to Oxygen. Fluid Phase Equilib., 19(3):175–200, 1985.
///     doi:10.1016/0378-3812(85)87016-3.
/// [2] Richard B. Stewart, Richard T. Jacobsen, and W. Wagner. Thermodynamic Properties of Oxygen
///     from the Triple Point to 300 K with Pressures to 80 MPa. J. Phys. Chem. Ref. Data,
///     20(5):917–1021, 1991. doi:10.1063/1.555897.
/// [3] E. W. Lemmon and R. T Jacobsen. Viscosity and Thermal Conductivity Equations for Nitrogen,
///     Oxygen, Argon, and Air. Int. J. Thermophys., 25(1):21–69, 2004.
///     doi:10.1023/B:IJOT.0000022327.04529.f3.
class Oxygen : public Helmholtz {
public:
    Oxygen();

private:
    [[nodiscard]] double alpha(double delta, double tau) const override;
    [[nodiscard]] double dalpha_ddelta(double delta, double tau) const override;
    [[nodiscard]] double dalpha_dtau(double delta, double tau) const override;
    [[nodiscard]] double d2alpha_ddelta2(double delta, double tau) const override;
    [[nodiscard]] double d2alpha_dtau2(double delta, double tau) const override;
    [[nodiscard]] double d2alpha_ddeltatau(double delta, double tau) const override;
    [[nodiscard]] double mu_from_rho_T(double rho, double T) const override;
    [[nodiscard]] double k_from_rho_T(double rho, double T) const override;

    IdealGasLead<double> lead;
    IdealGasLogTau<double> log_tau;
    IdealGasPlanckEinstein<double> pe;
    IdealEnthalpyEntropyOffset<double> ofst;
    ResidualPower<double> power_r;
    ResidualPowerExp<double, unsigned int> power_exp_r;

    CollisionIntegral<double> eta_0;
    ModifiedBatshinskiHildebrand<double> eta_r;
    Eta0AndPoly<double> lambda_0;
    PolynomialAndExponential<double> lambda_r;
};

} // namespace fprops
