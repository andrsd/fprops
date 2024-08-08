// SPDX-FileCopyrightText: 2022 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#pragma once

#include "Helmholtz.h"
#include "TransportModels.h"

namespace fprops {

/// Nitrogen (N2) fluid properties
///
/// References:
/// 1. Span,. Lemmon, Jacobsen, Wagner and Yokozeki, A reference equation of state for the
///    thermodynamic properties of nitrogen for temperatures from 63.151 to 1000 K and pressures to
///    2200 MPa, Journal of Physical and Chemical Reference Data, 29, 1361--1433 (2000)
/// 2. E. W. Lemmon and R. T Jacobsen. Viscosity and Thermal Conductivity Equations for Nitrogen,
///    Oxygen, Argon, and Air. Int. J. Thermophys., 25(1):21–69, 2004.
class Nitrogen : public Helmholtz {
public:
    Nitrogen();

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
    IdealGasPower<double> power_0;
    IdealGasPlanckEinsteinFunctionT<double> pefnt;
    ResidualPower<double> power_r;
    ResidualPowerExp<double, unsigned int> power_exp_r;
    ResidualGaussian<double> gauss;

    LennardJones<double> eta_0;
    ModifiedBatshinskiHildebrand<double> eta_r;
    Eta0AndPoly<double> lambda_0;
    ModifiedBatshinskiHildebrand<double> lambda_r;
};

} // namespace fprops
