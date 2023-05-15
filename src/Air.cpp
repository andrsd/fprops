#include "Air.h"
#include "Numerics.h"
#include <cmath>
#include <array>

namespace fprops {

const double two_thirds = 2. / 3.;

// Ref [1], Table 12
static std::array<double, 14> N = { 0,
                                    0.605719400e-7,
                                    -0.210274769e-4,
                                    -0.158860716e-3,
                                    -13.841928076,
                                    17.275266575,
                                    -0.195363420e-3,
                                    2.490888032,
                                    0.791309509,
                                    0.212236768,
                                    -0.197938904,
                                    25.36365,
                                    16.90741,
                                    87.31279 };

// Ref [1], Table 13
static std::array<double, 20> Nk = { 0,
                                     0.118160747229,
                                     0.713116392079,
                                     -0.161824192067e1,
                                     0.714140178971e-1,
                                     -0.865421396646e-1,
                                     0.134211176704,
                                     0.112626704218e-1,
                                     -0.420533228842e-1,
                                     0.349008431982e-1,
                                     0.164957183186e-3,
                                     -0.101365037912,
                                     -0.173813690970,
                                     -0.472103183731e-1,
                                     -0.122523554253e-1,
                                     -0.146629609713,
                                     -0.316055879821e-1,
                                     0.233594806142e-3,
                                     0.148287891978e-1,
                                     -0.938782884667e-2 };

const double a1 = 10.3753039487406;
const double a2 = 3.31112445645577;

// Ref [1], Table 13
static std::array<unsigned int, 20> ik = { 0, 1, 1, 1, 2, 3, 3, 4,  4, 4,
                                           6, 1, 3, 5, 6, 1, 3, 11, 1, 3 };

// Ref [1], Table 13
static std::array<double, 20> jk = { 0,    0,   0.33, 1.01, 0,    0,   0.15, 0,    0.2, 0.35,
                                     1.35, 1.6, 0.8,  0.95, 1.25, 3.6, 6,    3.25, 3.5, 15 };

// Ref [1], Table 13
static std::array<unsigned int, 20> lk = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 1, 1, 1, 1, 2, 2, 2, 3, 3 };

// Coefficients for viscosity
static std::array<double, 5> b_mu = { 0.431, -0.4623, 0.08406, 0.005341, -0.00331 };
static std::array<double, 5> N_mu = { 10.72, 1.122, 0.002019, -8.876, -0.02916 };
static std::array<double, 5> t_mu = { 0.2, 0.05, 2.4, 0.6, 3.6 };
static std::array<unsigned int, 5> d_mu = { 1, 4, 9, 1, 8 };
static std::array<unsigned int, 5> l_mu = { 0, 0, 0, 1, 1 };
static std::array<double, 5> gamma_mu = { 0.0, 0.0, 0.0, 1.0, 1.0 };

// Coefficients for thermal conductivity
static std::array<double, 6> N_k = { 8.743, 14.76, -16.62, 3.793, -6.142, -0.3778 };
static std::array<double, 6> t_k = { 0.1, 0.0, 0.5, 2.7, 0.3, 1.3 };
static std::array<unsigned int, 6> d_k = { 1, 2, 3, 7, 7, 11 };
static std::array<unsigned int, 6> l_k = { 0, 0, 2, 2, 2, 2 };
static std::array<double, 6> gamma_k = { 0.0, 0.0, 1.0, 1.0, 1.0, 1.0 };

Air::Air() : Helmholtz(8.314510, 28.96546e-3, 342.684564168, 132.5306) {}

double
Air::alpha(double delta, double tau) const
{
    // Ref [1], Eqn (24)
    double alpha_0 = std::log(delta);
    for (unsigned int i = 1; i <= 5; i++)
        alpha_0 += N[i] * std::pow(tau, i - 4);
    alpha_0 += N[6] * std::pow(tau, 1.5);
    alpha_0 += N[7] * std::log(tau);
    alpha_0 += N[8] * std::log(1 - std::exp(-N[11] * tau));
    alpha_0 += N[9] * std::log(1 - std::exp(-N[12] * tau));
    alpha_0 += N[10] * std::log(two_thirds + std::exp(N[13] * tau));
    alpha_0 += a1 + a2 * tau;

    // Ref [1], Eqn (26)
    double alpha_r = 0.0;
    for (unsigned int i = 1; i <= 10; i++)
        alpha_r += Nk[i] * std::pow(delta, ik[i]) * std::pow(tau, jk[i]);
    for (unsigned int i = 11; i <= 19; i++)
        alpha_r += Nk[i] * std::pow(delta, ik[i]) * std::pow(tau, jk[i]) *
                   std::exp(-std::pow(delta, lk[i]));

    return alpha_0 + alpha_r;
}

double
Air::dalpha_ddelta(double delta, double tau) const
{
    double dalpha_0 = 1.0 / delta;

    // Ref [1], Eqn (50)
    double dalpha_r = 0.0;
    for (unsigned int k = 1; k <= 10; k++)
        dalpha_r += ik[k] * Nk[k] * std::pow(delta, ik[k] - 1) * std::pow(tau, jk[k]);
    for (unsigned int k = 11; k <= 19; k++)
        dalpha_r += Nk[k] * std::pow(delta, ik[k] - 1) * std::pow(tau, jk[k]) *
                    std::exp(-std::pow(delta, lk[k])) * (ik[k] - lk[k] * std::pow(delta, lk[k]));

    return dalpha_0 + dalpha_r;
}

double
Air::dalpha_dtau(double delta, double tau) const
{
    // Ref [1], Eqn (48)
    double dalpha_0 = 0;
    for (unsigned i = 1; i <= 5; i++)
        dalpha_0 += (i - 4) * N[i] * std::pow(tau, i - 5);
    dalpha_0 += 1.5 * N[6] * std::pow(tau, 0.5);
    dalpha_0 += N[7] / tau;
    dalpha_0 += N[8] * N[11] / (std::exp(N[11] * tau) - 1);
    dalpha_0 += N[9] * N[12] / (std::exp(N[12] * tau) - 1);
    dalpha_0 += N[10] * N[13] / (two_thirds * std::exp(-N[13] * tau) + 1);
    dalpha_0 += a2;

    // Ref [1], Eqn (53)
    double dalpha_r = 0.0;
    for (unsigned int k = 1; k <= 10; k++)
        dalpha_r += jk[k] * Nk[k] * std::pow(delta, ik[k]) * std::pow(tau, jk[k] - 1);
    for (unsigned int k = 11; k <= 19; k++)
        dalpha_r += jk[k] * Nk[k] * std::pow(delta, ik[k]) * std::pow(tau, jk[k] - 1) *
                    std::exp(-std::pow(delta, lk[k]));

    return dalpha_0 + dalpha_r;
}

double
Air::d2alpha_ddelta2(double delta, double tau) const
{
    const double dalpha_0 = -1.0 / sqr(delta);

    // Ref [1], Eqn (51)
    double dalpha_r = 0.0;
    for (unsigned int k = 1; k <= 10; k++)
        dalpha_r += ik[k] * (ik[k] - 1) * Nk[k] * std::pow(delta, ik[k] - 2) * std::pow(tau, jk[k]);
    for (unsigned int k = 11; k <= 19; k++)
        dalpha_r += Nk[k] * std::pow(delta, ik[k] - 2) * std::pow(tau, jk[k]) *
                    std::exp(-std::pow(delta, lk[k])) *
                    ((ik[k] - lk[k] * std::pow(delta, lk[k])) *
                     (ik[k] - 1 - lk[k] * std::pow(delta, lk[k]) -
                      (lk[k] * lk[k] * std::pow(delta, lk[k]))));

    return dalpha_0 + dalpha_r;
}

double
Air::d2alpha_dtau2(double delta, double tau) const
{
    // Ref [1], Eqn (49)
    double dalpha_0 = 0.;
    for (unsigned int i = 1; i <= 5; i++)
        dalpha_0 += (i - 4) * (i - 5) * N[i] * std::pow(tau, i - 6);
    dalpha_0 += 0.75 * N[6] * std::pow(tau, -0.5);
    dalpha_0 -= N[7] / sqr(tau);
    dalpha_0 -= N[8] * sqr(N[11]) * std::exp(N[11] * tau) / sqr(std::exp(N[11] * tau) - 1);
    dalpha_0 -= N[9] * sqr(N[12]) * std::exp(N[12] * tau) / sqr(std::exp(N[12] * tau) - 1);
    dalpha_0 += two_thirds * N[10] * sqr(N[13]) * std::exp(-N[13] * tau) /
                sqr(two_thirds * std::exp(-N[13] * tau) + 1);

    // Ref [1], Eqn (54)
    double dalpha_r = 0.0;
    for (unsigned int k = 1; k <= 10; k++)
        dalpha_r += jk[k] * (jk[k] - 1) * Nk[k] * std::pow(delta, ik[k]) * std::pow(tau, jk[k] - 2);
    for (unsigned int k = 11; k <= 19; k++)
        dalpha_r += jk[k] * (jk[k] - 1) * Nk[k] * std::pow(delta, ik[k]) *
                    std::pow(tau, jk[k] - 2) * std::exp(-std::pow(delta, lk[k]));

    return dalpha_0 + dalpha_r;
}

double
Air::d2alpha_ddeltatau(double delta, double tau) const
{
    // Ref [1], Eqn (55)
    double dalpha_r = 0.0;
    for (unsigned int k = 1; k <= 10; k++)
        dalpha_r += ik[k] * jk[k] * Nk[k] * std::pow(delta, ik[k] - 1) * std::pow(tau, jk[k] - 1);
    for (unsigned int k = 11; k <= 19; k++)
        dalpha_r += jk[k] * Nk[k] * std::pow(delta, ik[k] - 1) * std::pow(tau, jk[k] - 1) *
                    std::exp(-std::pow(delta, lk[k])) * (ik[k] - lk[k] * std::pow(delta, lk[k]));

    return dalpha_r;
}

double
Air::eta0(double T) const
{
    double log_T_star = std::log(T / 103.3);
    double Omega_T_star = 0;
    for (unsigned int i = 0; i < b_mu.size(); i++)
        Omega_T_star += b_mu[i] * std::pow(log_T_star, i);
    Omega_T_star = std::exp(Omega_T_star);
    const double sigma = 0.360;
    return 0.0266958 * std::sqrt(1000.0 * this->M * T) / (sigma * sigma * Omega_T_star);
}

double
Air::mu_from_rho_T(double rho, double T) const
{
    double eta_0 = eta0(T);

    const double delta = rho / this->rho_c;
    const double tau = this->T_c / T;
    double eta_r = 0.0;
    for (unsigned int i = 0; i < N_mu.size(); i++)
        eta_r += N_mu[i] * std::pow(tau, t_mu[i]) * std::pow(delta, d_mu[i]) *
                 std::exp(-gamma_mu[i] * std::pow(delta, l_mu[i]));

    return (eta_0 + eta_r) * 1.0e-6;
}

double
Air::k_from_rho_T(double rho, double T) const
{
    const double delta = rho / this->rho_c;
    const double tau = this->T_c / T;

    double lambda_0 = 0.0;
    lambda_0 += 1.308 * eta0(T);
    lambda_0 += 1.405 * std::pow(tau, -1.1);
    lambda_0 += -1.036 * std::pow(tau, -0.3);

    double lambda_r = 0.0;
    for (unsigned int i = 0; i < N_k.size(); i++)
        lambda_r += N_k[i] * std::pow(tau, t_k[i]) * std::pow(delta, d_k[i]) *
                    std::exp(-gamma_k[i] * std::pow(delta, l_k[i]));

    return (lambda_0 + lambda_r) * 1.0e-3;
}

} // namespace fprops
