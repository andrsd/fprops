// SPDX-FileCopyrightText: 2022 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#pragma once

#include "Helmholtz.h"
#include "TransportModels.h"

namespace fprops {

/// Carbon dioxide fluid properties
///
/// References:
/// 1. R. Span and W. Wagner. A New Equation of State for Carbon Dioxide Covering the Fluid Region
///    from the Triple Point Temperature to 1100 K at Pressures up to 800 MPa. J. Phys. Chem. Ref.
///    Data, 25:1509–1596, 1996. doi:10.1063/1.555991.
/// 2. G. Scalabrin, P. Marchi, F. Finezzo, and R. Span. A Reference Multiparameter Thermal
///    Conductivity Equation for Carbon Dioxide with an Optimized Functional Form. J. Phys. Chem.
///    Ref. Data, 35(4):1549–1575, 2006. doi:10.1063/1.2213631.
/// 3. A. Fenghour, W.A. Wakeham, and V. Vesovic. The viscosity of carbon dioxide. J. Phys. Chem.
///    Ref. Data, 27(1):31–44, 1998. 5. doi:10.1063/1.556013.
class CarbonDioxide : public Helmholtz {
public:
    CarbonDioxide();

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
    IdealEnthalpyEntropyOffset<double> offset;

    ResidualPower<double> power_r;
    ResidualPowerExp<double, unsigned int> power_exp_r;
    ResidualGaussian<double> gauss;
    ResidualNonAnalytic<double> noan;

    LennardJones<double> eta_0;
    ModifiedBatshinskiHildebrand<double> eta_r;
    ModifiedBatshinskiHildebrand<double> lambda_r;
};

} // namespace fprops
