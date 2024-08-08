// SPDX-FileCopyrightText: 2022 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#include "fprops/carbon_dioxide.h"

namespace fprops {

static const double GAS_CONSTANT = 8.31451;
static const double MOLAR_MASS = 44.0098e-3;
static const double T_CRIT = 304.1282;
static const double RHO_MOLAR_CRIT = 10624.9063;
static const double RHO_CRIT = RHO_MOLAR_CRIT * MOLAR_MASS;

CarbonDioxide::CarbonDioxide() :
    Helmholtz(GAS_CONSTANT, MOLAR_MASS, RHO_CRIT, T_CRIT),
    lead(8.37304456, -3.70454304),
    log_tau(2.5),
    pe({ 1.99427042, 0.62105248, 0.41195293, 1.04028922, 0.08327678 },
       { 3.15163, 6.1119, 6.77708, 11.32384, 27.08792 }),
    offset(-14.4979156224319, 8.82013935801453),
    power_r({ 0.388568232032,
              2.93854759427,
              -5.5867188535,
              -0.767531995925,
              0.317290055804,
              0.548033158978,
              0.122794112203 },
            { 1, 1, 1, 1, 2, 2, 3 },
            { 0, 0.75, 1, 2, 0.75, 2, 0.75 }),
    power_exp_r(
        { 2.16589615432,      1.58417351097,    -0.231327054055,   0.0581169164314,
          -0.553691372054,    0.489466159094,   -0.0242757398435,  0.0624947905017,
          -0.121758602252,    -0.370556852701,  -0.0167758797004,  -0.11960736638,
          -0.0456193625088,   0.0356127892703,  -0.00744277271321, -0.00173957049024,
          -0.0218101212895,   0.0243321665592,  -0.0374401334235,  0.143387157569,
          -0.134919690833,    -0.0231512250535, 0.0123631254929,   0.00210583219729,
          -0.000339585190264, 0.00559936517716, -0.000303351180556 },
        { 1, 2, 4, 5, 5, 5, 6, 6, 6, 1, 1, 4, 4, 4, 7, 8, 2, 3, 3, 5, 5, 6, 7, 8, 10, 4, 8 },
        { 1.5, 1.5, 2.5, 0,  1.5, 2,  0,  1,  2,  3, 6, 3,  6, 8,
          6,   0,   7,   12, 16,  22, 24, 16, 24, 8, 2, 28, 14 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 6 }),
    gauss({ -213.654886883, 26641.5691493, -24027.2122046, -283.41603424, 212.472844002 },
          { 2, 2, 2, 3, 3 },
          { 1, 0, 1, 3, 3 },
          { 25, 25, 25, 15, 20 },
          { 1, 1, 1, 1, 1 },
          { 325, 300, 300, 275, 275 },
          { 1.16, 1.19, 1.19, 1.25, 1.22 }),
    noan({ -0.666422765408, 0.726086323499, 0.0550686686128 },
         { 3.5, 3.5, 3 },
         { 0.875, 0.925, 0.875 },
         { 0.3, 0.3, 0.3 },
         { 0.7, 0.7, 0.7 },
         { 0.3, 0.3, 1 },
         { 10, 10, 12.5 },
         { 275, 275, 275 }),

    eta_0(0.15178953643112786,
          MOLAR_MASS,
          251.196,
          1.,
          { 0.235156, -0.491266, 0.05211155, 0.05347906, -0.01537102 }),
    eta_r({ 1.9036541208525784,
            15.7384720473354,
            0.14207809578440784,
            0.0679058431241662,
            -0.030732988514867565 },
          { 0, 0, 3, 0, 1 },
          { 1, 2, 6, 8, 8 },
          { 0, 0, 0, 0, 0 },
          { 1, 1, 1, 1, 0 }),
    lambda_r({ 37.0597124660408,
               0.7696647124242399,
               7.5538113451464,
               -32.416436589336,
               78.894098855904,
               17.7830586854928,
               107.44756315137599,
               318.39746259479995,
               -0.82691726160072,
               2.0846013855224798e-02 },
             { 0.0, 0.0, -1.5, 0.0, -1.0, -1.5, -1.5, -1.5, -3.5, -5.5 },
             { 1.0, 5.0, 1.0, 1.0, 2.0, 0.0, 5.0, 9.0, 0.0, 0.0 },
             { 0, 0, 0, 2, 2, 2, 2, 2, 2, 2 },
             { 0, 0, 0, 5, 5, 5, 5, 5, 5, 5 })
{
}

double
CarbonDioxide::alpha(double delta, double tau) const
{
    // clang-format off
    return
        this->lead.alpha(delta, tau) +
        this->log_tau.alpha(delta, tau) +
        this->pe.alpha(delta, tau) +
        this->offset.alpha(delta, tau) +
        this->power_r.alpha(delta, tau) +
        this->power_exp_r.alpha(delta, tau) +
        this->gauss.alpha(delta, tau) +
        this->noan.alpha(delta, tau);
    // clang-format on
}

double
CarbonDioxide::dalpha_ddelta(double delta, double tau) const
{
    // clang-format off
    return
        this->lead.ddelta(delta, tau) +
        this->power_r.ddelta(delta, tau) +
        this->power_exp_r.ddelta(delta, tau) +
        this->gauss.ddelta(delta, tau) +
        this->noan.ddelta(delta, tau);
    // clang-format on
}

double
CarbonDioxide::dalpha_dtau(double delta, double tau) const
{
    // clang-format off
    return
        this->lead.dtau(delta, tau) +
        this->log_tau.dtau(delta, tau) +
        this->pe.dtau(delta, tau) +
        this->offset.dtau(delta, tau) +
        this->power_r.dtau(delta, tau) +
        this->power_exp_r.dtau(delta, tau) +
        this->gauss.dtau(delta, tau) +
        this->noan.dtau(delta, tau);
    // clang-format on
}

double
CarbonDioxide::d2alpha_ddelta2(double delta, double tau) const
{
    // clang-format off
    return
        this->lead.d2delta(delta, tau) +
        this->power_r.d2delta(delta, tau) +
        this->power_exp_r.d2delta(delta, tau) +
        this->gauss.d2delta(delta, tau) +
        this->noan.d2delta(delta, tau);
    // clang-format on
}

double
CarbonDioxide::d2alpha_dtau2(double delta, double tau) const
{
    // clang-format off
    return
        this->log_tau.d2tau(delta, tau) +
        this->pe.d2tau(delta, tau) +
        this->power_r.d2tau(delta, tau) +
        this->power_exp_r.d2tau(delta, tau) +
        this->gauss.d2tau(delta, tau) +
        this->noan.d2tau(delta, tau);
    // clang-format on
}

double
CarbonDioxide::d2alpha_ddeltatau(double delta, double tau) const
{
    // clang-format off
    return
        this->power_r.d2deltatau(delta, tau) +
        this->power_exp_r.d2deltatau(delta, tau) +
        this->gauss.d2deltatau(delta, tau) +
        this->noan.d2deltatau(delta, tau);
    // clang-format on
}

double
CarbonDioxide::mu_from_rho_T(double rho, double T) const
{
    const double d = delta(rho);
    const double t = tau(T);

    double eta = this->eta_0.value(T) + this->eta_r.value(d, t);
    // [Pa-s]
    return eta * 1.0e-6;
}

double
CarbonDioxide::k_from_rho_T(double rho, double T) const
{
    const double d = delta(rho);
    const double t = tau(T);

    double lambda = this->lambda_r.value(d, t);
    // [W/(m-K)]
    return lambda * 1.0e-3;
}

} // namespace fprops
