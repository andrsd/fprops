// SPDX-FileCopyrightText: 2024 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#include "fprops/oxygen.h"

namespace fprops {

static const double GAS_CONSTANT = 8.31434;
static const double MOLAR_MASS = 31.9988e-3;
static const double T_CRIT = 154.581;
static const double RHO_MOLAR_CRIT = 13630.0;
static const double RHO_CRIT = RHO_MOLAR_CRIT * MOLAR_MASS;

Oxygen::Oxygen() :
    Helmholtz(GAS_CONSTANT, MOLAR_MASS, RHO_CRIT, T_CRIT),
    lead(0, 0),
    log_tau(2.51808732),
    pe({ 1.02323928, 0.784357918, 0.00337183363, -0.0170864084, 0.0463751562 },
       { 14.5316979447668,
         72.8419165356674,
         7.7710849975094,
         0.446425786480874,
         34.4677188658373 }),
    ofst(-14.716836666461498, -0.011083985429237124),
    power_r({ 0.3983768749,
              -1.846157454,
              0.4183473197,
              0.02370620711,
              0.09771730573,
              0.03017891294,
              0.02273353212,
              0.01357254086,
              -0.04052698943,
              0.0005454628515,
              0.0005113182277,
              2.953466883e-07,
              -8.687645072e-05 },
            { 1, 1, 1, 2, 2, 2, 3, 3, 3, 6, 7, 7, 8 },
            { 0, 1.5, 2.5, -0.5, 1.5, 2, 0, 1, 2.5, 0, 2, 5, 2 }),
    power_exp_r({ -0.2127082589,
                  0.08735941958,
                  0.127550919,
                  -0.09067701064,
                  -0.03540084206,
                  -0.03623278059,
                  0.0132769929,
                  -0.0003254111865,
                  -0.008313582932,
                  0.002124570559,
                  -0.0008325206232,
                  -2.626173276e-05,
                  0.002599581482,
                  0.009984649663,
                  0.002199923153,
                  -0.02591350486,
                  -0.1259630848,
                  0.1478355637,
                  -0.01011251078 },
                { 1, 1, 2, 2, 3, 3, 5, 6, 7, 8, 10, 2, 3, 3, 4, 4, 5, 5, 5 },
                { 5, 6, 3.5, 5.5, 3, 7, 6, 8.5, 4, 6.5, 5.5, 22, 11, 18, 11, 23, 17, 18, 23 },
                { 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4 }),
    eta_0(2.66958e-08,
          0.0319988,
          118.5,
          3.428e-10,
          { 0.431, -0.4623, 0.08406, 0.005341, -0.00331 },
          { 0, 1, 2, 3, 4 }),
    eta_r({ 1.767e-05, 4.042e-07, 1.077e-10, 3.51e-07, -1.367e-05 },
          { 1, 5, 12, 8, 1 },
          { 0.05, 0.0, 2.1, 0.0, 0.5 },
          { 0, 0, 0, -1, -1 },
          { 0, 0, 0, 1, 2 },
          { 0 },
          { 1 },
          { 0 },
          { 1 },
          { 0 },
          { 1 },
          { 0 }),
    lambda_0({ 0.001036, 0.006283, -0.004262 }, { 0, -0.9, -0.6 }),
    lambda_r({ 0.01531, 0.008898, -0.0007336, 0.006728, -0.004374, -0.0004747 },
             { 0.0, 0.0, 0.3, 4.3, 0.5, 1.8 },
             { 1, 3, 4, 5, 7, 10 },
             { 0, 0, 0, 1, 1, 1 },
             { 0, 0, 0, 2, 2, 2 })
{
}

double
Oxygen::alpha(double delta, double tau) const
{
    // clang-format off
    return
        this->lead.alpha(delta, tau) +
        this->log_tau.alpha(delta, tau) +
        this->pe.alpha(delta, tau) +
        this->ofst.alpha(delta, tau) +
        this->power_r.alpha(delta, tau) +
        this->power_exp_r.alpha(delta, tau);
    // clang-format on
}

double
Oxygen::dalpha_ddelta(double delta, double tau) const
{
    // clang-format off
    return
        this->lead.ddelta(delta, tau) +
        this->power_r.ddelta(delta, tau) +
        this->power_exp_r.ddelta(delta, tau);
    // clang-format on
}

double
Oxygen::dalpha_dtau(double delta, double tau) const
{
    // clang-format off
    return
        this->lead.dtau(delta, tau) +
        this->log_tau.dtau(delta, tau) +
        this->pe.dtau(delta, tau) +
        this->ofst.dtau(delta, tau) +
        this->power_r.dtau(delta, tau) +
        this->power_exp_r.dtau(delta, tau);
    // clang-format on
}

double
Oxygen::d2alpha_ddelta2(double delta, double tau) const
{
    // clang-format off
    return
        this->lead.d2delta(delta, tau) +
        this->power_r.d2delta(delta, tau) +
        this->power_exp_r.d2delta(delta, tau);
    // clang-format on
}

double
Oxygen::d2alpha_dtau2(double delta, double tau) const
{
    // clang-format off
    return
        this->log_tau.d2tau(delta, tau) +
        this->pe.d2tau(delta, tau) +
        this->power_r.d2tau(delta, tau) +
        this->power_exp_r.d2tau(delta, tau);
    // clang-format on
}

double
Oxygen::d2alpha_ddeltatau(double delta, double tau) const
{
    // clang-format off
    return
        this->power_r.d2deltatau(delta, tau) +
        this->power_exp_r.d2deltatau(delta, tau);
    // clang-format on
}

double
Oxygen::mu_from_rho_T(double rho, double T) const
{
    auto d = delta(rho);
    auto t = tau(T);
    // clang-format off
    return
        this->eta_0.value(T) +
        this->eta_r.value(d, t);
    // clang-format on
}

double
Oxygen::k_from_rho_T(double rho, double T) const
{
    auto d = delta(rho);
    auto t = tau(T);
    double eta0 = this->eta_0.value(T);
    // clang-format off
    return
        this->lambda_0.value(eta0, t) +
        this->lambda_r.value(t, d);
    // clang-format on
}

} // namespace fprops
