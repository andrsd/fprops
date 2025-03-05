// SPDX-FileCopyrightText: 2022 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#pragma once

#include "fprops/helmholtz.h"
#include "fprops/transport_models.h"

namespace fprops {

/// Air fluid properties
///
/// References:
/// 1. Eric W. Lemmon, Richard T. Jacobsen, Steven G. Penoncello, and Daniel G. Friend.
///    Thermodynamic Properties of Air and Mixtures of Nitrogen, Argon, and Oxygen from 60 to
///    2000 K at Pressures to 2000 MPa. J. Phys. Chem. Ref. Data, 29(3):331–385, 2000.
/// 2. E. W. Lemmon and R. T Jacobsen. Viscosity and Thermal Conductivity Equations for Nitrogen,
///    Oxygen, Argon, and Air. Int. J. Thermophys., 25(1):21–69, 2004.
class Air : public Helmholtz<Air> {
public:
    Air();

private:
    [[nodiscard]] double alpha(double delta, double tau) const;
    [[nodiscard]] double dalpha_ddelta(double delta, double tau) const;
    [[nodiscard]] double dalpha_dtau(double delta, double tau) const;
    [[nodiscard]] double d2alpha_ddelta2(double delta, double tau) const;
    [[nodiscard]] double d2alpha_dtau2(double delta, double tau) const;
    [[nodiscard]] double d2alpha_ddeltatau(double delta, double tau) const;
    [[nodiscard]] double mu_from_rho_T(double rho, double T) const;
    [[nodiscard]] double k_from_rho_T(double rho, double T) const;

    IdealGasLead<double> lead;
    IdealGasPower<double> power_0;
    IdealGasLogTau<double> log_tau;
    IdealGasPlanckEinstein<double> pe;
    IdealPlanckEinsteinGeneralized<double> pegen;
    IdealEnthalpyEntropyOffset<double> offset;
    ResidualPower<double> power_r;
    ResidualPowerExp<double, unsigned int> power_exp_r;

    CollisionIntegral<double> eta_0;
    ModifiedBatshinskiHildebrand<double> eta_r;
    Eta0AndPoly<double> lambda_0;
    PolynomialAndExponential<double> lambda_r;

    template <typename FLUID>
    friend class Helmholtz;
};

} // namespace fprops
