// SPDX-FileCopyrightText: 2024 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#pragma once

#include "fprops/helmholtz.h"
#include "fprops/transport_models.h"

namespace fprops {

/// Methane fluid properties
///
/// Reference:
///
class Methane : public Helmholtz<Methane> {
public:
    Methane();

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
    IdealGasPlanckEinsteinFunctionT<double> pefnt;
    IdealEnthalpyEntropyOffset<double> offset;
    ResidualPower<double> power_r;
    ResidualPowerExp<double, unsigned int> power_exp_r;
    ResidualGaussian<double> gauss;

    PowersOfTemperature<double> eta_0;
    FrictionTheory<double> eta_f;

    template <typename FLUID>
    friend class Helmholtz;
};

} // namespace fprops
