// SPDX-FileCopyrightText: 2024 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#pragma once

#include "fprops/helmholtz.h"
#include "fprops/transport_models.h"

namespace fprops {

/// Ammonia fluid properties
///
class Ammonia : public Helmholtz<Ammonia> {
public:
    Ammonia();

private:
    [[nodiscard]] double alpha(double delta, double tau) const;
    [[nodiscard]] double dalpha_ddelta(double delta, double tau) const;
    [[nodiscard]] double dalpha_dtau(double delta, double tau) const;
    [[nodiscard]] double d2alpha_ddelta2(double delta, double tau) const;
    [[nodiscard]] double d2alpha_dtau2(double delta, double tau) const;
    [[nodiscard]] double d2alpha_ddeltatau(double delta, double tau) const;
    [[nodiscard]] double mu_from_rho_T(double rho, double T) const;
    [[nodiscard]] double k_from_rho_T(double rho, double T) const;

    double lambda_crit(double t, double d) const;

    IdealGasLead<double> lead;
    IdealGasLogTau<double> log_tau;
    IdealGasPlanckEinstein<double> pe;
    ResidualPower<double> power_r;
    ResidualPowerExp<double, unsigned int> power_exp_r;
    ResidualGaussian<double> gauss;
    ResidualGaoB<double> gaob;

    CollisionIntegral<double> eta_0;
    RainwaterFriend<double> rf;
    ModifiedBatshinskiHildebrand<double> eta_ho;
    PolynomialRatio<double> lambda_0;
    Polynomial<double> lambda_r;

    template <typename FLUID>
    friend class Helmholtz;
};

} // namespace fprops
