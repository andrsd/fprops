// SPDX-FileCopyrightText: 2022 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#include "fprops/air.h"

namespace fprops {

static const double GAS_CONSTANT = 8.314510;
static const double MOLAR_MASS = 28.96546e-3;
static const double T_CRIT = 132.5306;
static const double RHO_MOLAR_CRIT = 11830.8;
static const double RHO_CRIT = RHO_MOLAR_CRIT * MOLAR_MASS;

Air::Air() :
    Helmholtz(GAS_CONSTANT, MOLAR_MASS, RHO_CRIT, T_CRIT),
    lead(0, 0),
    power_0({ 6.057194e-08,
              -2.10274769e-05,
              -0.000158860716,
              -13.841928076,
              17.275266575,
              -0.00019536342 },
            { -3, -2, -1, 0, 1, 1.5 }),
    log_tau(2.490888032),
    pe({ 0.791309509, 0.212236768 }, { 25.36365, 16.90741 }),
    pegen({ -0.197938904 }, { 87.31279 }, { 0.6666666666666666 }, { 1 }),
    offset(10.3753039487406, 3.31112445645577),
    power_r({ 0.118160747229,
              0.713116392079,
              -1.61824192067,
              0.0714140178971,
              -0.0865421396646,
              0.134211176704,
              0.0112626704218,
              -0.0420533228842,
              0.0349008431982,
              0.000164957183186 },
            { 1, 1, 1, 2, 3, 3, 4, 4, 4, 6 },
            { 0, 0.33, 1.01, 0, 0, 0.15, 0, 0.2, 0.35, 1.35 }),
    power_exp_r({ -0.101365037912,
                  -0.17381369097,
                  -0.0472103183731,
                  -0.0122523554253,
                  -0.146629609713,
                  -0.0316055879821,
                  0.000233594806142,
                  0.0148287891978,
                  -0.00938782884667 },
                { 1, 3, 5, 6, 1, 3, 11, 1, 3 },
                { 1.6, 0.8, 0.95, 1.25, 3.6, 6, 3.25, 3.5, 15 },
                { 1, 1, 1, 1, 2, 2, 2, 3, 3 }),
    eta_0(2.66958e-08,
          0.0289586,
          103.3,
          3.6e-10,
          { 0.431, -0.4623, 0.08406, 0.005341, -0.00331 },
          { 0, 1, 2, 3, 4 }),
    eta_r({ 1.072e-05, 1.122e-06, 2.019e-09, -8.876e-06, -2.916e-08 },
          { 1, 4, 9, 1, 8 },
          { 0.2, 0.05, 2.4, 0.6, 3.6 },
          { 0, 0, 0, -1, -1 },
          { 0, 0, 0, 1, 1 },
          { 0 },
          { 1 },
          { 0 },
          { 1 },
          { 0 },
          { 1 },
          { 0 }),
    lambda_0({ 0.001308, 0.001405, -0.001036 }, { 0, -1.1, -0.3 }),
    lambda_r({ 0.008743, 0.01476, -0.01662, 0.003793, -0.006142, -0.0003778 },
             { 0.1, 0.0, 0.5, 2.7, 0.3, 1.3 },
             { 1, 2, 3, 7, 7, 11 },
             { 0.0, 0.0, 1.0, 1.0, 1.0, 1.0 },
             { 0, 0, 2, 2, 2, 2 })
{
}

double
Air::alpha(double delta, double tau) const
{
    // clang-format off
    return
        this->lead.alpha(delta, tau) +
        this->power_0.alpha(delta, tau) +
        this->log_tau.alpha(delta, tau) +
        this->pe.alpha(delta, tau) +
        this->pegen.alpha(delta, tau) +
        this->offset.alpha(delta, tau) +
        this->power_r.alpha(delta, tau) +
        this->power_exp_r.alpha(delta, tau);
    // clang-format on
}

double
Air::dalpha_ddelta(double delta, double tau) const
{
    // clang-format off
    return
        this->lead.ddelta(delta, tau) +
        this->power_r.ddelta(delta, tau) +
        this->power_exp_r.ddelta(delta, tau);
    // clang-format on
}

double
Air::dalpha_dtau(double delta, double tau) const
{
    // clang-format off
    return
        this->lead.dtau(delta, tau) +
        this->power_0.dtau(delta, tau) +
        this->log_tau.dtau(delta, tau) +
        this->pe.dtau(delta, tau) +
        this->pegen.dtau(delta, tau) +
        this->offset.dtau(delta, tau) +
        this->power_r.dtau(delta, tau) +
        this->power_exp_r.dtau(delta, tau);
    // clang-format on
}

double
Air::d2alpha_ddelta2(double delta, double tau) const
{
    // clang-format off
    return
        this->lead.d2delta(delta, tau) +
        this->power_r.d2delta(delta, tau) +
        this->power_exp_r.d2delta(delta, tau);
    // clang-format on
}

double
Air::d2alpha_dtau2(double delta, double tau) const
{
    // clang-format off
    return
        this->power_0.d2tau(delta, tau) +
        this->log_tau.d2tau(delta, tau) +
        this->pe.d2tau(delta, tau) +
        this->pegen.d2tau(delta, tau) +
        this->power_r.d2tau(delta, tau) +
        this->power_exp_r.d2tau(delta, tau);
    // clang-format on
}

double
Air::d2alpha_ddeltatau(double delta, double tau) const
{
    // clang-format off
    return
        this->power_r.d2deltatau(delta, tau) +
        this->power_exp_r.d2deltatau(delta, tau);
    // clang-format on
}

double
Air::mu_from_rho_T(double rho, double T) const
{
    auto d = delta(rho);
    auto t = tau(T);
    double eta = this->eta_0.value(T) + this->eta_r.value(d, t);
    return eta;
}

double
Air::k_from_rho_T(double rho, double T) const
{
    const double d = delta(rho);
    const double t = tau(T);
    double eta0 = this->eta_0.value(T);
    // FIXME: add the critical part
    return this->lambda_0.value(eta0, t) + this->lambda_r.value(t, d);
}

} // namespace fprops
