// SPDX-FileCopyrightText: 2022 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#include "fprops/nitrogen.h"

namespace fprops {

static const double MOLAR_MASS = 28.01348e-3;
static const double T_CRIT = 126.192;

Nitrogen::Nitrogen() :
    Helmholtz(8.314510, MOLAR_MASS, 313.299958972, T_CRIT),
    lead(-12.76952708, -0.00784163),
    log_tau(2.5),
    power_0({ -0.0001934819, -1.247742e-05, 6.678326e-08 }, { -1, -2, -3 }),
    pefnt(T_CRIT, { 1.012941 }, { 3364.011 }),
    power_r({ 0.924803575275,
              -0.492448489428,
              0.661883336938,
              -1.92902649201,
              -0.0622469309629,
              0.349943957581 },
            { 1, 1, 2, 2, 3, 3 },
            { 0.25, 0.875, 0.5, 0.875, 0.375, 0.75 }),
    power_exp_r({ 0.564857472498,    -1.61720005987,     -0.481395031883,   0.421150636384,
                  -0.0161962230825,  0.172100994165,     0.00735448924933,  0.0168077305479,
                  -0.00107626664179, -0.0137318088513,   0.000635466899859, 0.00304432279419,
                  -0.0435762336045,  -0.0723174889316,   0.0389644315272,   -0.021220136391,
                  0.00408822981509,  -5.51990017984e-05, -0.0462016716479,  -0.00300311716011,
                  0.0368825891208,   -0.0025585684622,   0.00896915264558,  -0.0044151337035,
                  0.00133722924858,  0.000264832491957 },
                { 1, 1, 1, 3, 3, 4, 6, 6, 7, 7, 8, 8, 1, 2, 3, 4, 5, 8, 4, 5, 5, 8, 3, 5, 6, 9 },
                { 0.5, 0.75, 2, 1.25, 3.5, 1,  0.5, 3,  0,  2.75, 0.75, 2.5, 4,
                  6,   6,    3, 3,    6,   16, 11,  15, 12, 12,   7,    4,   16 },
                { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 }),
    gauss({ 19.6688194015, -20.911560073, 0.0167788306989, 2627.67566274 },
          { 1, 1, 3, 2 },
          { 0, 1, 2, 3 },
          { 20, 20, 15, 25 },
          { 1, 1, 1, 1 },
          { 325, 325, 300, 275 },
          { 1.16, 1.16, 1.13, 1.25 }),

    eta_0(2.66958e-08,
          0.02801348,
          98.94,
          3.656e-10,
          { 0.431, -0.4623, 0.08406, 0.005341, -0.00331 },
          { 0, 1, 2, 3, 4 }),
    eta_r({ 1.072e-05, 3.989e-08, 1.208e-09, -7.402e-06, 4.62e-06 },
          { 2, 10, 12, 2, 1 },
          { 0.1, 0.25, 3.2, 0.9, 0.3 },
          { 0, -1, -1, -1, -1 },
          { 0, 1, 1, 2, 3 },
          { 0 },
          { 1 },
          { 0 },
          { 1 },
          { 0 },
          { 1 },
          { 0 }),
    lambda_0({ 0.001511, 0.002117, -0.003332 }, { 0, -1, -0.7 }),
    lambda_r({ 0.008862, 0.03111, -0.07313, 0.02003, -0.0007096, 0.0002672 },
             { 0.0, 0.03, 0.2, 0.8, 0.6, 1.9 },
             { 1, 2, 3, 4, 8, 10 },
             { 0.0, 0.0, 1.0, 1.0, 1.0, 1.0 },
             { 0, 0, 1, 2, 2, 2 })
{
}

double
Nitrogen::alpha(double delta, double tau) const
{
    // clang-format off
    return
        this->lead.alpha(delta, tau) +
        this->log_tau.alpha(delta, tau) +
        this->power_0.alpha(delta, tau) +
        this->pefnt.alpha(delta, tau) +

        this->power_r.alpha(delta, tau) +
        this->power_exp_r.alpha(delta, tau) +
        this->gauss.alpha(delta, tau);
    // clang-format on
}

double
Nitrogen::dalpha_ddelta(double delta, double tau) const
{
    // clang-format off
    return
        this->lead.ddelta(delta, tau) +

        this->power_r.ddelta(delta, tau) +
        this->power_exp_r.ddelta(delta, tau) +
        this->gauss.ddelta(delta, tau);
    // clang-format on
}

double
Nitrogen::dalpha_dtau(double delta, double tau) const
{
    // clang-format off
    return
        this->lead.dtau(delta, tau) +
        this->log_tau.dtau(delta, tau) +
        this->power_0.dtau(delta, tau) +
        this->pefnt.dtau(delta, tau) +

        this->power_r.dtau(delta, tau) +
        this->power_exp_r.dtau(delta, tau) +
        this->gauss.dtau(delta, tau);
    // clang-format on
}

double
Nitrogen::d2alpha_ddelta2(double delta, double tau) const
{
    // clang-format off
    return
        this->lead.d2delta(delta, tau) +

        this->power_r.d2delta(delta, tau) +
        this->power_exp_r.d2delta(delta, tau) +
        this->gauss.d2delta(delta, tau);
    // clang-format on
}

double
Nitrogen::d2alpha_dtau2(double delta, double tau) const
{
    // clang-format off
    return
        this->log_tau.d2tau(delta, tau) +
        this->power_0.d2tau(delta, tau) +
        this->pefnt.d2tau(delta, tau) +

        this->power_r.d2tau(delta, tau) +
        this->power_exp_r.d2tau(delta, tau) +
        this->gauss.d2tau(delta, tau);
    // clang-format on
}

double
Nitrogen::d2alpha_ddeltatau(double delta, double tau) const
{
    // clang-format off
    return
        this->power_r.d2deltatau(delta, tau) +
        this->power_exp_r.d2deltatau(delta, tau) +
        this->gauss.d2deltatau(delta, tau);
    // clang-format on
}

double
Nitrogen::mu_from_rho_T(double rho, double T) const
{
    const double d = delta(rho);
    const double t = tau(T);
    return this->eta_0.value(T) + this->eta_r.value(d, t);
}

double
Nitrogen::k_from_rho_T(double rho, double T) const
{
    const double d = delta(rho);
    const double t = tau(T);

    double eta0 = this->eta_0.value(T);
    // FIXME: add the critical part
    return this->lambda_0.value(eta0, t) + this->lambda_r.value(t, d);
}

} // namespace fprops
