// SPDX-FileCopyrightText: 2022 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#pragma once

#include "fprops/helmholtz.h"
#include "fprops/transport_models.h"

namespace fprops {

/// Helium (He) fluid properties
///
/// References:
/// 1. D. O. Ortiz-Vega, K. R. Hall, J. C. Holste, V. D. Arp, A. H. Harvey, and E. W. Lemmon.
///    Equation of state for Helium-4. Unpublished.
/// 2. B.A. Hands and V.D. Arp. A Correlation of Thermal Conductivity Data for Helium. Cryogenics,
///    21(12):697–703, 1981. doi:10.1016/0011-2275(81)90211-3.
/// 3. V.D. Arp, R.D. McCarty, and D.G Friend. Thermophysical Properties of Helium-4 from 0.8 to
///    1500 K with Pressures to 2000 MPa - NIST Technical Note 1334 (revised). Technical Report,
///    NIST, 1998.
class Helium : public Helmholtz<Helium> {
public:
    Helium();

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
    IdealEnthalpyEntropyOffset<double> offset;
    ResidualPower<double> power_r;
    ResidualPowerExp<double, unsigned int> power_exp_r;
    ResidualGaussian<double> gauss;

    template <typename FLUID>
    friend class Helmholtz;
};

} // namespace fprops
