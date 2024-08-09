// SPDX-FileCopyrightText: 2024 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#include "fprops/ammonia.h"
#include <iomanip>
#include <iostream>

namespace fprops {

static const double GAS_CONSTANT = 8.3144598;
static const double MOLAR_MASS = 17.03052e-3;
static const double T_CRIT = 405.56;
static const double RHO_MOLAR_CRIT = 13696;
static const double RHO_CRIT = RHO_MOLAR_CRIT * MOLAR_MASS;
static const double EPSILON_OVER_K = 386;
static const double SIGMA_ETA = 2.957e-10;

Ammonia::Ammonia() :
    Helmholtz(GAS_CONSTANT, MOLAR_MASS, RHO_CRIT, T_CRIT),
    lead(-6.59406093943886, 5.60101151987913),
    log_tau(3.),
    pe({ 2.224, 3.148, 0.9579 }, { 4.0585856593352405, 9.776605187888352, 17.829667620080876 }),
    power_r({ 0.006132232, 1.7395866, -2.2261792, -0.30127553, 0.08967023 },
            { 4, 1, 1, 2, 3 },
            { 1.0, 0.382, 1.0, 1.0, 0.677 }),
    power_exp_r({ -0.076387037, -0.84063963, -0.27026327 },
                { 3, 2, 3 },
                { 2.915, 3.51, 1.063 },
                { 2, 2, 1 }),
    gauss({ 6.212578,
            -5.7844357,
            2.4817542,
            -2.3739168,
            0.01493697,
            -3.7749264,
            0.0006254348,
            -1.7359e-05,
            -0.13462033,
            0.07749072839 },
          { 1, 1, 1, 2, 2, 1, 3, 3, 1, 1 },
          { 0.655, 1.3, 3.1, 1.4395, 1.623, 0.643, 1.13, 4.5, 1.0, 4.0 },
          { 0.42776, 0.6424, 0.8175, 0.7995, 0.91, 0.3574, 1.21, 4.14, 22.56, 22.68 },
          { -0.0726, -0.1274, 0.7527, 0.57, 2.2, -0.243, 2.96, 3.02, 0.9574, 0.9576 },
          { 1.708, 1.4865, 2.0915, 2.43, 0.488, 1.1, 0.85, 1.14, 945.64, 993.85 },
          { 1.036, 1.2777, 1.083, 1.2906, 0.928, 0.934, 0.919, 1.852, 1.05897, 1.05277 }),
    gaob({ -1.6909858, 0.93739074 },
         { 1, 1 },
         { 4.3315, 4.015 },
         { -2.8452, -2.8342 },
         { 0.3696, 0.2962 },
         { 1.108, 1.313 },
         { 0.4478, 0.44689 },
         { 1.244, 0.6826 }),
    //
    eta_0(2.1357e-06,
          0.01703026,
          EPSILON_OVER_K,
          SIGMA_ETA,
          { 4.9931822, -0.61122364, 0.0, 0.18535124, -0.11160946 },
          { 0, 1, 2, 3, 4 }),
    rf(EPSILON_OVER_K,
       SIGMA_ETA,
       { -1.7999496,
         46.692621,
         -534.60794,
         3360.4074,
         -13019.164,
         33414.23,
         -58711.743,
         71426.686,
         -59834.012,
         33652.741,
         -12027.35,
         2434.8205,
         -208.07957 },
       { -0.0, -0.5, -1.0, -1.5, -2.0, -2.5, -3.0, -3.5, -4.0, -4.5, -5.0, -5.5, -6.0 }),
    eta_ho({ 4.005040600989671e-06,
             -1.4107915123955129e-05,
             3.4760743039321816e-05,
             4.631310990138071e-06,
             -3.937374461785061e-06,
             -1.200075068367531e-05,
             1.9284977991745303e-06 },
           { 3, 3, 2, 4, 4, 2, 4 },
           { 0, 1, 2, 2, 3, 4, 4 },
           { 0, 0, 0, 0, 0, 0, 0 },
           { 1, 1, 1, 1, 1, 1, 1 },
           { 0.0 },
           { 1 },
           { 0 },
           { 1 },
           { 0 },
           { 1 },
           { 0 }),
    lambda_0({ 0.03589, -0.000175, 4.551e-07, 1.685e-10, -4.828e-13 },
             { 0, 1, 2, 3, 4 },
             { 1.0 },
             { 0 }),
    lambda_r({ 0.03808645, 0.06647986, -0.0300295, 0.00998779 }, { 0, 0, 0, 0 }, { 1, 2, 3, 4 })
{
}

double
Ammonia::alpha(double delta, double tau) const
{
    // clang-format off
    return
        this->lead.alpha(delta, tau) +
        this->log_tau.alpha(delta, tau) +
        this->pe.alpha(delta, tau) +
        this->power_r.alpha(delta, tau) +
        this->power_exp_r.alpha(delta, tau) +
        this->gauss.alpha(delta, tau) +
        this->gaob.alpha(delta, tau);
    // clang-format on
}

double
Ammonia::dalpha_ddelta(double delta, double tau) const
{
    // clang-format off
    return
        this->lead.ddelta(delta, tau) +
        this->power_r.ddelta(delta, tau) +
        this->power_exp_r.ddelta(delta, tau) +
        this->gauss.ddelta(delta, tau) +
        this->gaob.ddelta(delta, tau);
    // clang-format on
}

double
Ammonia::dalpha_dtau(double delta, double tau) const
{
    // clang-format off
    return
        this->lead.dtau(delta, tau) +
        this->log_tau.dtau(delta, tau) +
        this->pe.dtau(delta, tau) +
        this->power_r.dtau(delta, tau) +
        this->power_exp_r.dtau(delta, tau) +
        this->gauss.dtau(delta, tau) +
        this->gaob.dtau(delta, tau);
    // clang-format on
}

double
Ammonia::d2alpha_ddelta2(double delta, double tau) const
{
    // clang-format off
    return
        this->lead.d2delta(delta, tau) +
        this->power_r.d2delta(delta, tau) +
        this->power_exp_r.d2delta(delta, tau) +
        this->gauss.d2delta(delta, tau) +
        this->gaob.d2delta(delta, tau);
    // clang-format on
}

double
Ammonia::d2alpha_dtau2(double delta, double tau) const
{
    // clang-format off
    return
        this->log_tau.d2tau(delta, tau) +
        this->pe.d2tau(delta, tau) +
        this->power_r.d2tau(delta, tau) +
        this->power_exp_r.d2tau(delta, tau) +
        this->gauss.d2tau(delta, tau) +
        this->gaob.d2tau(delta, tau);
    // clang-format on
}

double
Ammonia::d2alpha_ddeltatau(double delta, double tau) const
{
    // clang-format off
    return
        this->power_r.d2deltatau(delta, tau) +
        this->power_exp_r.d2deltatau(delta, tau) +
        this->gauss.d2deltatau(delta, tau) +
        this->gaob.d2deltatau(delta, tau);
    // clang-format on
}

double
Ammonia::mu_from_rho_T(double rho, double T) const
{
    const double d = delta(rho);
    const double t = tau(T);

    auto eta_dilute = this->eta_0.value(T);
    auto B_eta_initial = this->rf.value(T);
    auto rho_molar = rho / MOLAR_MASS;
    auto initial_density = eta_dilute * B_eta_initial * rho_molar;
    // clang-format off
    return
        eta_dilute +
        initial_density +
        this->eta_ho.value(d, t);
    // clang-format on
}

double
Ammonia::k_from_rho_T(double rho, double T) const
{
    double t = 1. / T;
    double d = rho / 235.;
    // clang-format off
    return
        this->lambda_0.value(T) +
        this->lambda_r.value(t, d) +
        lambda_crit(rho, T);
    // clang-format on
}

double
Ammonia::lambda_crit(double rho, double T) const
{
    // From "Thermal Conductivity of Ammonia in a Large Temperature and Pressure Range Including the
    // Critical Region" by R. Tufeu, D.Y. Ivanov, Y. Garrabos, B. Le Neindre, Bereicht der
    // Bunsengesellschaft Phys. Chem. 88 (1984) 422-427

    double Tc = 405.4;
    double rhoc = 235;
    double LAMBDA = 1.2;
    double nu = 0.63;
    double gamma = 1.24;
    double DELTA = 0.50;
    double zeta_0_plus = 1.34e-10;
    double a_zeta = 1;
    double GAMMA_0_plus = 0.423e-8;
    double pi = 3.141592654;
    double k_B = 1.3806504e-23;

    auto t = std::abs((T - Tc) / Tc);
    auto a_chi = a_zeta / 0.7;
    auto eta_B = (2.60 + 1.6 * t) * 1e-5;
    auto dPdT = (2.18 - 0.12 / std::exp(17.8 * t)) * 1e5; // [Pa-K]
    auto X_T = 0.61 * rhoc + 16.5 * std::log(t);
    // Along the critical isochore (only a function of temperature) (Eq. 9)
    auto DELTA_lambda_i =
        LAMBDA * (k_B * T * T) /
        (6 * pi * eta_B * (zeta_0_plus * pow(t, -nu) * (1 + a_zeta * pow(t, DELTA)))) * dPdT *
        dPdT * GAMMA_0_plus * math::pow(t, -gamma) * (1 + a_chi * pow(t, DELTA));
    auto DELTA_lambda_id = DELTA_lambda_i * std::exp(-36 * t * t);
    double DELTA_lambda;
    if (rho < 0.6 * rhoc) {
        DELTA_lambda = DELTA_lambda_id * (X_T * X_T) /
                       (X_T * X_T + math::pow<2>(0.6 * rhoc - 0.96 * rhoc)) * math::pow<2>(rho) /
                       math::pow<2>(0.6 * rhoc);
    }
    else {
        DELTA_lambda =
            DELTA_lambda_id * (X_T * X_T) / (X_T * X_T + math::pow<2>(rho - 0.96 * rhoc));
    }
    return DELTA_lambda;
}

} // namespace fprops
