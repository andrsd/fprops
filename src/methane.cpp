// SPDX-FileCopyrightText: 2024 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#include "fprops/methane.h"

namespace fprops {

static const double GAS_CONSTANT = 8.314510;
static const double MOLAR_MASS = 16.0428e-3;
static const double T_CRIT = 190.564;
static const double RHO_MOLAR_CRIT = 10139.128;
static const double RHO_CRIT = RHO_MOLAR_CRIT * MOLAR_MASS;

Methane::Methane() :
    Helmholtz(GAS_CONSTANT, MOLAR_MASS, RHO_CRIT, T_CRIT),
    lead(9.91243972, -6.33270087),
    log_tau(3.0016),
    pefnt(T_CRIT, { 0.008449, 4.6942, 3.4865, 1.6572, 1.4115 }, { 648, 1957, 3895, 5705, 15080 }),
    offset(-12.8829893867948, 9.22344625310864),
    power_r({ 0.04367901028,
              0.6709236199,
              -1.765577859,
              0.8582330241,
              -1.206513052,
              0.512046722,
              -0.0004000010791,
              -0.01247842423,
              0.03100269701,
              0.001754748522,
              -3.171921605e-06,
              -2.24034684e-06,
              2.947056156e-07 },
            { 1, 1, 1, 2, 2, 2, 2, 3, 4, 4, 8, 9, 10 },
            { -0.5, 0.5, 1, 0.5, 1, 1.5, 4.5, 0, 1, 3, 1, 3, 3 }),
    power_exp_r({ 0.1830487909,   0.1511883679,   -0.4289363877,  0.06894002446, -0.01408313996,
                  -0.0306305483,  -0.02969906708, -0.01932040831, -0.1105739959, 0.09952548995,
                  0.008548437825, -0.06150555662, -0.04291792423, -0.0181320729, 0.0344590476,
                  -0.00238591945, -0.01159094939, 0.06641693602,  -0.0237154959, -0.03961624905,
                  -0.01387292044, 0.03389489599,  -0.002927378753 },
                { 1, 1, 1, 2, 4, 5, 6, 1, 2, 3, 4, 4, 3, 5, 5, 8, 2, 3, 4, 4, 4, 5, 6 },
                { 0, 1, 2, 0, 0, 2, 2, 5, 5, 5, 2, 4, 12, 8, 10, 10, 10, 14, 12, 18, 22, 18, 14 },
                { 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4 }),
    gauss({ 9.324799946e-05, -6.287171518, 12.71069467, -6.423953466 },
          { 2, 0, 0, 0 },
          { 2, 0, 1, 2 },
          { 20, 40, 40, 40 },
          { 1, 1, 1, 1 },
          { 200, 250, 250, 250 },
          { 1.07, 1.11, 1.11, 1.11 }),
    //
    eta_0({ 2.60536e-06, -1.85247e-05, 2.34216e-05, 0.0 }, { 0, 0.25, 0.5, 0.75 }),
    eta_f(1,
          1,
          { 3.49668e-08, -1.73176e-08, 0.0 },
          { -3.12118e-08, 1.99422e-10, 0 },
          1,
          { 5.98858e-08, -4.91143e-08, 0.0 },
          1,
          { -8.52992e-13, -3.58009e-13, 0.0 },
          3,
          {},
          { 1.60099e-11, 8.50221e-13, 0.0 },
          3,
          { -3.55631e-10, 2.80326e-10, 0.0 },
          3,
          {},
          0,
          {},
          0)
{
}

double
Methane::alpha(double delta, double tau) const
{
    // clang-format off
    return
        this->lead.alpha(delta, tau) +
        this->log_tau.alpha(delta, tau) +
        this->pefnt.alpha(delta, tau) +
        this->offset.alpha(delta, tau) +
        this->power_r.alpha(delta, tau) +
        this->power_exp_r.alpha(delta, tau) +
        this->gauss.alpha(delta, tau);
    // clang-format on
}

double
Methane::dalpha_ddelta(double delta, double tau) const
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
Methane::dalpha_dtau(double delta, double tau) const
{
    // clang-format off
    return
        this->lead.dtau(delta, tau) +
        this->log_tau.dtau(delta, tau) +
        this->pefnt.dtau(delta, tau) +
        this->offset.dtau(delta, tau) +
        this->power_r.dtau(delta, tau) +
        this->power_exp_r.dtau(delta, tau) +
        this->gauss.dtau(delta, tau);
    // clang-format on
}

double
Methane::d2alpha_ddelta2(double delta, double tau) const
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
Methane::d2alpha_dtau2(double delta, double tau) const
{
    // clang-format off
    return
        this->log_tau.d2tau(delta, tau) +
        this->pefnt.d2tau(delta, tau) +
        this->power_r.d2tau(delta, tau) +
        this->power_exp_r.d2tau(delta, tau) +
        this->gauss.d2tau(delta, tau);
    // clang-format on
}

double
Methane::d2alpha_ddeltatau(double delta, double tau) const
{
    // clang-format off
    return
        this->power_r.d2deltatau(delta, tau) +
        this->power_exp_r.d2deltatau(delta, tau) +
        this->gauss.d2deltatau(delta, tau);
    // clang-format on
}

double
Methane::mu_from_rho_T(double rho, double T) const
{
    auto Tr = T / T_CRIT;

    auto d = delta(rho);
    auto t = tau(T);
    auto da_dd = dalpha_ddelta(d, t);
    auto p = pressure(rho, T, d, da_dd);
    // [bar]; 1e5 for conversion from Pa -> bar
    auto p_bar = p / 1e5;
    auto rho_molar = rho / MOLAR_MASS;
    auto dp_dT = GAS_CONSTANT * rho_molar * d * da_dd;
    auto p_r_bar = T * dp_dT / 1e5;
    auto p_a_bar = p_bar - p_r_bar;
    auto p_id_bar = rho_molar * GAS_CONSTANT * T / 1e5;
    // clang-format off
    return
        this->eta_0.value(Tr) +
        this->eta_f.value(t, p_a_bar, p_r_bar, p_id_bar);
    // clang-format on
}

double
Methane::k_from_rho_T(double rho, double T) const
{
    const double d = delta(rho);
    const double t = tau(T);

    // NOTE: slightly different rho_crit and T_crit used in computation of
    // `eta`s
    double delta_eta = rho / 10139.0;
    double tau_eta = 190.55 / T;

    // Viscosity formulation from Friend, JPCRD, 1989
    // Dilute
    double C[] = { 0,
                   -3.0328138281,
                   16.918880086,
                   -37.189364917,
                   41.288861858,
                   -24.615921140,
                   8.9488430959,
                   -1.8739245042,
                   0.20966101390,
                   -9.6570437074e-3 };
    double OMEGA22_sum = 0;
    double t_eta = T / 174.0;
    for (int i = 1; i <= 9; ++i)
        OMEGA22_sum += C[i] * math::pow(t_eta, (i - 1.0) / 3.0 - 1.0);
    double eta_dilute = 10.50 * sqrt(t_eta) * OMEGA22_sum;
    // double re[] = { 0, 1, 1, 2, 2, 2, 3, 3, 4, 4, 1, 1 };
    // double se[] = { 0, 0, 1, 0, 1, 1.5, 0, 2, 0, 1, 0, 1 };
    // double ge[] = { 0,           0.41250137, -0.14390912, 0.10366993,  0.40287464,  -0.24903524,
    //                 -0.12953131, 0.06575776, 0.02566628,  -0.03716526, -0.38798341, 0.03533815 };
    // double sum1 = 0;
    // for (int i = 1; i <= 9; ++i)
    //     sum1 += ge[i] * math::pow(delta_eta, re[i]) * math::pow(tau_eta, se[i]);
    // double sum2 = 0;
    // for (int i = 10; i <= 11; ++i)
    //     sum2 += ge[i] * math::pow(delta_eta, re[i]) * math::pow(tau_eta, se[i]);
    //
    // double eta_residual = 12.149 * sum1 / (1 + sum2);
    // double eta = eta_residual + eta_dilute;

    // Dilute
    auto d2alpha0_dtau2 = this->log_tau.d2tau(d, t) + this->pefnt.d2tau(d, t);

    double f_int = 1.458850 - 0.4377162 / t_eta;
    auto lambda_dilute = 0.51828 * eta_dilute *
                         (3.75 - f_int * (math::pow<2>(tau(T)) * d2alpha0_dtau2 + 1.5)); // [mW/m/K]
    // Residual
    double rl[] = { 0, 1, 3, 4, 4, 5, 5, 2 };
    double sl[] = { 0, 0, 0, 0, 1, 0, 1, 0 };
    double jl[] = { 0,           2.4149207,  0.55166331,   -0.52837734,
                    0.073809553, 0.24465507, -0.047613626, 1.5554612 };
    double sum = 0;
    for (int i = 1; i <= 6; ++i) {
        sum += jl[i] * math::pow(delta_eta, rl[i]) * math::pow(tau_eta, sl[i]);
    }
    double delta_sigma_star = 1.0; // Looks like a typo in Friend - should be 1 instead of 11
    // FIXME: uncomment
    // if (HEOS.T() < HEOS.T_critical() && HEOS.rhomolar() < HEOS.rhomolar_critical()) {
    //     delta_sigma_star = HEOS.saturation_ancillary(iDmolar, 1, iT, HEOS.T()) /
    //                        HEOS.keyed_output(CoolProp::irhomolar_critical);
    // }
    auto lambda_residual =
        6.29638 * (sum + jl[7] * math::pow<2>(delta_eta) / delta_sigma_star); // [mW/m/K]

    // FIXME: Ignoring critical region
    // // Critical region
    // double Tstar = 1 - 1 / tau;
    // double rhostar = 1 - delta;
    // double F_T = 2.646, F_rho = 2.678, F_A = -0.637;
    // double F = std::exp(-F_T * std::sqrt(std::abs(Tstar)) - F_rho * math::pow<2>(rhostar) - F_A *
    // rhostar);
    // double CHI_T_star;
    // if (std::abs(Tstar) < 0.03) {
    //     if (std::abs(rhostar) < 1e-16) {
    //         // Equation 26
    //         const double LAMBDA = 0.0801, gamma = 1.190;
    //         CHI_T_star = LAMBDA * math::pow(std::abs(Tstar), -gamma);
    //     }
    //     else if (std::abs(rhostar) < 0.03) {
    //         // Equation 23
    //         const double beta = 0.355, W = -1.401, S = -6.098, E = 0.287, a = 3.352, b = 0.732,
    //                      R = 0.535, Q = 0.1133;
    //         double OMEGA = W * Tstar * math::pow(std::abs(rhostar), -1 / beta);
    //         double theta = 1;
    //         if (Tstar < -math::pow(std::abs(rhostar), -1 / beta) / S) {
    //             theta = 1 + E * math::pow(1 + S * Tstar * math::pow(std::abs(rhostar), -1 /
    //             beta), 2 * beta);
    //         }
    //         CHI_T_star =
    //             Q * math::pow(std::abs(rhostar), -a) * math::pow(theta, b) / (theta + OMEGA *
    //             (theta + R));
    //     }
    //     else {
    //         // Equation 19a
    //         CHI_T_star = 0.28631 * delta * tau /
    //                      (1 + 2 * delta * HEOS.dalphar_dDelta() +
    //                       math::pow<2>(delta) * HEOS.d2alphar_dDelta2());
    //     }
    // }
    // else {
    //     // Equation 19a
    //     CHI_T_star =
    //         0.28631 * delta * tau /
    //         (1 + 2 * delta * HEOS.dalphar_dDelta() + math::pow<2>(delta) *
    //         HEOS.d2alphar_dDelta2());
    // }

    // auto lambda_critical = 91.855 / (eta * math::pow<2>(tau)) *
    //                        math::pow<2>(1 + delta * HEOS.dalphar_dDelta() -
    //                                     delta * tau * HEOS.d2alphar_dDelta_dTau()) *
    //                        math::pow(CHI_T_star, 0.4681) * F; //[mW/m/K]
    // auto lambda = (lambda_dilute + lambda_residual + lambda_critical) * 0.001;

    return (lambda_dilute + lambda_residual) * 0.001;
}

} // namespace fprops
