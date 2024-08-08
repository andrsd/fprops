// SPDX-FileCopyrightText: 2022 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#include "fprops/helium.h"

namespace fprops {

static const double GAS_CONSTANT = 8.3144598;
static const double MOLAR_MASS = 4.002602e-3;
static const double T_CRIT = 5.1953;
static const double RHO_MOLAR = 18130.0;
static const double RHO_CRIT = RHO_MOLAR * MOLAR_MASS;

Helium::Helium() :
    Helmholtz(GAS_CONSTANT, MOLAR_MASS, RHO_CRIT, T_CRIT),
    lead(0.1871304489697973, 0.4848903984696551),
    log_tau(1.5),
    offset(-0.0278608575394308, -0.00835828466224539),
    power_r({ 0.015559018, 3.0638932, -4.2420844, 0.054418088, -0.18971904, 0.087856262 },
            { 4, 1, 1, 2, 2, 3 },
            { 1.0, 0.425, 0.63, 0.69, 1.83, 0.575 }),
    power_exp_r({ 2.2833566, -0.53331595, -0.53296502, 0.99444915, -0.30078896, -1.6432563 },
                { 1, 1, 3, 2, 2, 1 },
                { 0.925, 1.585, 1.69, 1.51, 2.9, 0.8 },
                { 1, 2, 2, 1, 2, 1 }),
    gauss({ 0.8029102,
            0.026838669,
            0.04687678,
            -0.14832766,
            0.03016211,
            -0.019986041,
            0.14283514,
            0.007418269,
            -0.22989793,
            0.79224829,
            -0.049386338 },
          { 2, 1, 2, 1, 1, 3, 2, 2, 3, 2, 2 },
          { 1.26, 3.51, 2.785, 1.0, 4.22, 0.83, 1.575, 3.447, 0.73, 1.634, 6.13 },
          { 1.5497, 9.245, 4.76323, 6.3826, 8.7023, 0.255, 0.3523, 0.1492, 0.05, 0.1668, 42.2358 },
          { 0.596, 0.3423, 0.761, 0.9747, 0.5868, 0.5627, 2.5346, 3.6763, 4.5245, 5.039, 0.959 },
          { 0.2471,
            0.0983,
            0.1556,
            2.6782,
            2.7077,
            0.6621,
            0.1775,
            0.4821,
            0.3069,
            0.1758,
            1357.6577 },
          { 3.15, 2.54505, 1.2513, 1.9416, 0.5984, 2.2282, 1.606, 3.815, 1.61958, 0.6407, 1.076 })
{
}

double
Helium::alpha(double delta, double tau) const
{
    // clang-format off
    return
        this->lead.alpha(delta, tau) +
        this->log_tau.alpha(delta, tau) +
        this->offset.alpha(delta, tau) +

        this->power_r.alpha(delta, tau) +
        this->power_exp_r.alpha(delta, tau) +
        this->gauss.alpha(delta, tau);
    // clang-format on
}

double
Helium::dalpha_ddelta(double delta, double tau) const
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
Helium::dalpha_dtau(double delta, double tau) const
{
    // clang-format off
    return
        this->lead.dtau(delta, tau) +
        this->log_tau.dtau(delta, tau) +
        this->offset.dtau(delta, tau) +

        this->power_r.dtau(delta, tau) +
        this->power_exp_r.dtau(delta, tau) +
        this->gauss.dtau(delta, tau);
    // clang-format on
}

double
Helium::d2alpha_ddelta2(double delta, double tau) const
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
Helium::d2alpha_dtau2(double delta, double tau) const
{
    // clang-format off
    return
        this->log_tau.d2tau(delta, tau) +

        this->power_r.d2tau(delta, tau) +
        this->power_exp_r.d2tau(delta, tau) +
        this->gauss.d2tau(delta, tau);
    // clang-format on
}

double
Helium::d2alpha_ddeltatau(double delta, double tau) const
{
    // clang-format off
    return
        this->power_r.d2deltatau(delta, tau) +
        this->power_exp_r.d2deltatau(delta, tau) +
        this->gauss.d2deltatau(delta, tau);
    // clang-format on
}

double
Helium::mu_from_rho_T(double rho, double T) const
{
    // This is taken from CoolProp

    // Correlation uses density in g/cm^3, ours is kg/m^3
    rho = rho / 1000.0;

    double x;
    if (T <= 300)
        x = std::log(T);
    else
        x = std::log(300.0);

    // Evaluate the terms B, C, D
    double B = -47.5295259 / x + 87.6799309 - 42.0741589 * x + 8.33128289 * math::pow<2>(x) -
               0.589252385 * math::pow<3>(x);
    double C = 547.309267 / x - 904.870586 + 431.404928 * x - 81.4504854 * math::pow<2>(x) +
               5.37008433 * math::pow<3>(x);
    double D = -1684.39324 / x + 3331.08630 - 1632.19172 * x + 308.804413 * math::pow<2>(x) -
               20.2936367 * math::pow<3>(x);
    double eta_0_slash = -0.135311743 / x + 1.00347841 + 1.20654649 * x -
                         0.149564551 * math::pow<2>(x) + 0.012520841 * math::pow<3>(x);
    double eta_E_slash = rho * B + rho * rho * C + rho * rho * rho * D;

    double eta;
    if (T <= 100) {
        double ln_eta = eta_0_slash + eta_E_slash;
        eta = std::exp(ln_eta);
    }
    else {
        double ln_eta = eta_0_slash + eta_E_slash;
        double eta_0 = 196 * std::pow(T, 0.71938) * std::exp(12.451 / T - 295.67 / T / T - 4.1249);
        eta = std::exp(ln_eta) + eta_0 - std::exp(eta_0_slash);
    }
    // Correlation yields viscosity in micro g/(cm-s); divide by 10 to get micro Pa-s
    eta = eta / 10.;
    // [Pa-s]
    return eta * 1.0e-6;
}

double
Helium::k_from_rho_T(double rho, double T) const
{
    // This is taken from CoolProp
    double rhoc = 68.0;
    double summer = 3.739232544 / T - 2.620316969e1 / T / T + 5.982252246e1 / T / T / T -
                    4.926397634e1 / T / T / T / T;
    double lambda_0 = 2.7870034e-3 * std::pow(T, 7.034007057e-1) * std::exp(summer);
    double c[] = { 1.862970530e-4,  -7.275964435e-7, -1.427549651e-4, 3.290833592e-5,
                   -5.213335363e-8, 4.492659933e-8,  -5.924416513e-9, 7.087321137e-6,
                   -6.013335678e-6, 8.067145814e-7,  3.995125013e-7 };
    // Equation 17
    double lambda_e =
        (c[0] + c[1] * T + c[2] * std::pow(T, 1 / 3.0) + c[3] * std::pow(T, 2.0 / 3.0)) * rho +
        (c[4] + c[5] * std::pow(T, 1.0 / 3.0) + c[6] * std::pow(T, 2.0 / 3.0)) * rho * rho * rho +
        (c[7] + c[8] * std::pow(T, 1.0 / 3.0) + c[9] * std::pow(T, 2.0 / 3.0) + c[10] / T) * rho *
            rho * std::log(rho / rhoc);

    // Critical component
    double lambda_c = 0.0;

    if ((3.5 < T) && (T < 12)) {
        double x0 = 0.392, E1 = 2.8461, E2 = 0.27156, beta = 0.3554, gamma = 1.1743, delta = 4.304,
               rhoc_crit = 69.158, Tc = 5.18992, pc = 2.2746e5;

        double DeltaT = std::abs(1 - T / Tc);
        double DeltaRho = std::abs(1 - rho / rhoc_crit);
        double eta = mu_from_rho_T(rho, T); // [Pa-s]

        // K_T = HEOS.isothermal_compressibility();
        // K_T = 1.0 / _rhomolar * first_partial_deriv(iDmolar, iP, iT);
        // K_T = 1.0 / drho_dpT()
        double K_T = 0; // FIXME
        double K_Tbar;
        double W = std::pow(DeltaT / 0.2, 2) + std::pow(DeltaRho / 0.25, 2);
        if (W > 1) {
            K_Tbar = K_T;
        }
        else {
            double x = std::pow(DeltaT / DeltaRho, 1 / beta);
            double h = E1 * (1 + x / x0) *
                       std::pow(1 + E2 * std::pow(1 + x / x0, 2 / beta), (gamma - 1) / (2 * beta));
            double dhdx =
                E1 *
                (E2 * std::pow((x + x0) / x0, 2 / beta) * (gamma - 1) *
                     std::pow(E2 * std::pow((x + x0) / x0, 2 / beta) + 1,
                              (1.0 / 2.0) * (gamma - 1) / beta) +
                 std::pow(beta, 2) * std::pow(E2 * pow((x + x0) / x0, 2 / beta) + 1,
                                              (1.0 / 2.0) * (2 * beta + gamma - 1) / beta)) /
                (std::pow(beta, 2) * x0 * (E2 * std::pow((x + x0) / x0, 2 / beta) + 1));
            // Right-hand-side of Equation 9
            double RHS = std::pow(DeltaRho, delta - 1) * (delta * h - x / beta * dhdx);
            double K_Tprime = 1 / (RHS * std::pow(rho / rhoc_crit, 2) * pc);
            K_Tbar = W * K_T + (1 - W) * K_Tprime;
        }

        // double dpdT = HEOS.first_partial_deriv(CoolProp::iP, CoolProp::iT, CoolProp::iDmolar);
        double dpdT = 0; // FIXME

        // 3.4685233d-17 and 3.726229668d0 are "magical" coefficients that are present in the
        // REFPROP source to yield the right values. Not clear why these values are needed. Also,
        // the form of the critical term in REFPROP does not agree with Hands paper. EL and MH from
        // NIST are not sure where these coefficients come from.
        lambda_c = 3.4685233e-17 * 3.726229668 * std::sqrt(K_Tbar) * std::pow(T, 2) / rho / eta *
                   std::pow(dpdT, 2) *
                   std::exp(-18.66 * std::pow(DeltaT, 2) - 4.25 * std::pow(DeltaRho, 4));
    }
    return lambda_0 + lambda_e + lambda_c;
}

} // namespace fprops
