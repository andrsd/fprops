// SPDX-FileCopyrightText: 2026 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#include "fprops/r134a.h"

namespace fprops {

static const double GAS_CONSTANT = 8.314471;
static const double MOLAR_MASS = 0.102032;
static const double T_CRIT = 374.21;
static const double RHO_MOLAR = 5017.053;
static const double RHO_CRIT = RHO_MOLAR * MOLAR_MASS;

R134a::R134a() :
    Helmholtz(GAS_CONSTANT, MOLAR_MASS, RHO_CRIT, T_CRIT),
    lead(-1.019535, 9.047135),
    log_tau(-1.629789),
    power({ -9.723916, -3.92717 }, { -0.5, -0.75 }),
    power_r({ 0.05586817,
              0.498223,
              0.02458698,
              0.0008570145,
              0.0004788584,
              -1.800808,
              0.2671641,
              -0.04781652 },
            { 2, 1, 3, 6, 6, 1, 1, 2 },
            { -0.5, 0, 0, 0, 1.5, 1.5, 2, 2 }),
    power_exp_r(
        // n
        { 0.01423987,
          0.3324062,
          -0.007485907,
          0.0001017263,
          -0.5184567,
          -0.08692288,
          0.2057144,
          -0.005000457,
          0.0004603262,
          -0.003497836,
          0.006995038,
          -0.01452184,
          -0.0001285458 },
        // d
        { 5, 2, 2, 4, 1, 4, 1, 2, 4, 1, 5, 3, 10 },
        // t
        { 1, 3, 5, 1, 5, 5, 6, 10, 10, 10, 18, 22, 50 },
        // l
        { 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4 }),
    //
    eta_0(2.1357e-08,
          0.102031,
          299.363,
          4.68932e-10,
          { 0.355404, -0.464337, 0.0257353 },
          { 0, 1, 2 }),
    eta_r({ -2.06900719e-05,
            3.56029549e-07,
            2.11101816e-06,
            1.39601415e-05,
            -4.5643502e-06,
            -3.51593275e-06 },
          { 1, 2, 2, 2, 2, 3 },
          { 0, 6, 2, 0.5, -2, 0 },
          { 0, 0, 0, 0, 0, 0 },
          { 1, 1, 1, 1, 1, 1 },
          { 0.00021476332 },
          { 0 },
          { 0 },
          { 3.163695636 },
          { 0 },
          { 1, -0.0890173375, 0.100035295 },
          { 0, -1, -2 }),
    //
    lambda_0({ -0.0105248, 8.00982e-05 }, { 0, 1 }, { 1.0 }, { 0 }),
    lambda_r({ 0.0037740609300000003, 0.010534223865, -0.002952794565, 0.00128672592 },
             { 0, 0, 0, 0 },
             { 1, 2, 3, 4 }),
    T_reducing(1.),
    rhomass_reducing(515.249968352)
{
}

double
R134a::alpha(double delta, double tau) const
{
    // clang-format off
    return
        this->lead.alpha(delta, tau) +
        this->log_tau.alpha(delta, tau) +

        this->power.alpha(delta, tau) +
        this->power_r.alpha(delta, tau) +
        this->power_exp_r.alpha(delta, tau)
    ;
    // clang-format on
}

double
R134a::dalpha_ddelta(double delta, double tau) const
{
    // clang-format off
    return
        this->lead.ddelta(delta, tau) +

        this->power_r.ddelta(delta, tau) +
        this->power_exp_r.ddelta(delta, tau)
    ;
    // clang-format on
}

double
R134a::dalpha_dtau(double delta, double tau) const
{
    // clang-format off
    return
        this->lead.dtau(delta, tau) +
        this->log_tau.dtau(delta, tau) +

        this->power.dtau(delta, tau) +
        this->power_r.dtau(delta, tau) +
        this->power_exp_r.dtau(delta, tau)
    ;
    // clang-format on
}

double
R134a::d2alpha_ddelta2(double delta, double tau) const
{
    // clang-format off
    return
        this->lead.d2delta(delta, tau) +

        this->power_r.d2delta(delta, tau) +
        this->power_exp_r.d2delta(delta, tau)
    ;
    // clang-format on
}

double
R134a::d2alpha_dtau2(double delta, double tau) const
{
    // clang-format off
    return
        this->log_tau.d2tau(delta, tau) +

        this->power.d2tau(delta, tau) +
        this->power_r.d2tau(delta, tau) +
        this->power_exp_r.d2tau(delta, tau)
    ;
    // clang-format on
}

double
R134a::d2alpha_ddeltatau(double delta, double tau) const
{
    // clang-format off
    return
        this->power_r.d2deltatau(delta, tau) +
        this->power_exp_r.d2deltatau(delta, tau)
    ;
    // clang-format on
}

double
R134a::mu_from_rho_T(double rho, double T) const
{
    auto d = delta(rho);
    auto t = tau(T);
    double eta = this->eta_0.value(T) + this->eta_r.value(d, t);
    return eta;
}

double
R134a::k_from_rho_T(double rho, double T) const
{
    double t = this->T_reducing / T;
    double Tr = T / this->T_reducing;
    const double d = rho / this->rhomass_reducing;
    // clang-format off
    return
        this->lambda_0.value(Tr) +
        this->lambda_r.value(t, d)
    ;
    // clang-format on
}

} // namespace fprops
