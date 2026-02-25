// SPDX-FileCopyrightText: 2026 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#pragma once

#include "fprops/helmholtz.h"
#include "fprops/transport_models.h"

namespace fprops {

class R134a : public Helmholtz<R134a> {
public:
    R134a();

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
    IdealGasLogTau<double> log_tau;
    IdealGasPower<double> power;
    ResidualPower<double> power_r;
    ResidualPowerExp<double, unsigned int> power_exp_r;
    // viscosity
    CollisionIntegral<double> eta_0;
    ModifiedBatshinskiHildebrand<double> eta_r;
    // conductivity
    PolynomialRatio<double> lambda_0;
    Polynomial<double> lambda_r;

    const double T_reducing;
    const double rhomass_reducing;

    template <typename FLUID>
    friend class Helmholtz;
};

} // namespace fprops
