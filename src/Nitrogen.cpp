#include "Nitrogen.h"
#include "Numerics.h"
#include <cmath>

namespace fprops {

#define len(a) (sizeof(a) / sizeof(a[0]))

#define R 8.314510
#define M 28.01348e-3
#define rho_c 313.299958972
#define T_c 126.192

// Coefficients for ideal gas component of the Helmholtz free energy
static const double a[] = { 2.5,          -12.76952708, -0.00784163, -1.934819e-4,
                            -1.247742e-5, 6.678326e-8,  1.012941,    26.65788 };

// Coefficients for residual component of the Helmholtz free energy
static const double N1[] = { 0.924803575275,    -0.492448489428,    0.661883336938,
                             -0.192902649201e1, -0.622469309629e-1, 0.349943957581 };
static const unsigned int i1[] = { 1, 1, 2, 2, 3, 3 };
static const double j1[] = { 0.25, 0.875, 0.5, 0.875, 0.375, 0.75 };

static const double N2[] = { 0.564857472498,     -0.161720005987e1,  -0.481395031883,
                             0.421150636384,     -0.161962230825e-1, 0.172100994165,
                             0.735448924933e-2,  0.168077305479e-1,  -0.107626664179e-2,
                             -0.137318088513e-1, 0.635466899859e-3,  0.304432279419e-2,
                             -0.435762336045e-1, -0.723174889316e-1, 0.389644315272e-1,
                             -0.21220136391e-1,  0.4808822981509e-2, -0.551990017984e-4,
                             -0.462016716479e-1, -0.300311716011e-2, 0.368825891208e-1,
                             -0.25585684622e-2,  0.896915264558e-2,  -0.44151337035e-2,
                             0.133722924858e-2,  0.264832491957e-3 };
static const unsigned int i2[] = { 1, 1, 1, 3, 3, 4, 6, 6, 7, 7, 8, 8, 1,
                                   2, 3, 4, 5, 8, 4, 5, 5, 8, 3, 5, 6, 9 };
static const double j2[] = { 0.5,  0.75, 2.0,  1.25, 3.5,  1.0, 0.5, 3.0, 0.0,
                             2.75, 0.75, 2.5,  4.0,  6.0,  6.0, 3.0, 3.0, 6.0,
                             16.0, 11.0, 15.0, 12.0, 12.0, 7.0, 4.0, 16.0 };
static const unsigned int l2[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2,
                                   2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 };

static const double N3[] = { 0.196688194015e2,
                             -0.20911560073e2,
                             0.167788306989e-1,
                             0.262767566274e4 };
static const unsigned int i3[] = { 1, 1, 3, 2 };
static const unsigned int j3[] = { 0, 1, 2, 3 };
static const unsigned int l3[] = { 2, 2, 2, 2 };
static const double phi3[] = { 20.0, 20.0, 15.0, 25.0 };
static const double beta3[] = { 325.0, 325.0, 300.0, 275.0 };
static const double gamma3[] = { 1.16, 1.16, 1.13, 1.25 };

/// Coefficients for viscosity
static const double b_mu[] = { 0.431, -0.4623, 0.08406, 0.005341, -0.00331 };
static const double N_mu[] = { 10.72, 0.03989, 0.001208, -7.402, 4.62 };
static const double t_mu[] = { 0.1, 0.25, 3.2, 0.9, 0.3 };
static const double d_mu[] = { 2, 10, 12, 2, 1 };
static const double l_mu[] = { 0, 1, 1, 2, 3 };
static const double gamma_mu[] = { 0.0, 1.0, 1.0, 1.0, 1.0 };

/// Coefficients for thermal conductivity
static const double N_k[] = { 8.862, 31.11, -73.13, 20.03, -0.7096, 0.2672 };
static const double t_k[] = { 0.0, 0.03, 0.2, 0.8, 0.6, 1.9 };
static const unsigned int d_k[] = { 1, 2, 3, 4, 8, 10 };
static const unsigned int l_k[] = { 0, 0, 1, 2, 2, 2 };
static const double gamma_k[] = { 0.0, 0.0, 1.0, 1.0, 1.0 };

Nitrogen::Nitrogen() : Helmholtz(R, M, rho_c, T_c) {}

double
Nitrogen::alpha(double delta, double tau)
{
    // Ideal gas component of the Helmholtz free energy
    double alpha0 = std::log(delta) + a[0] * std::log(tau) + a[1] + a[2] * tau + a[3] / tau +
                    a[4] / sqr(tau) + a[5] / cb(tau) + a[6] * std::log(1.0 - std::exp(-a[7] * tau));

    // Residual component of the Helmholtz free energy
    double alphar = 0.0;
    for (unsigned int i = 0; i < len(N1); i++)
        alphar += N1[i] * std::pow(delta, i1[i]) * std::pow(tau, j1[i]);
    for (unsigned int i = 0; i < len(N2); i++)
        alphar += N2[i] * std::pow(delta, i2[i]) * std::pow(tau, j2[i]) *
                  std::exp(-std::pow(delta, l2[i]));
    for (unsigned int i = 0; i < len(N3); i++)
        alphar += N3[i] * std::pow(delta, i3[i]) * std::pow(tau, j3[i]) *
                  std::exp(-phi3[i] * sqr(delta - 1.0) - beta3[i] * sqr(tau - gamma3[i]));

    return alpha0 + alphar;
}

double
Nitrogen::dalpha_ddelta(double delta, double tau)
{
    // Ideal gas component of the Helmholtz free energy
    double dalpha0 = 1.0 / delta;

    // Residual component of the Helmholtz free energy
    double dalphar = 0.0;
    for (unsigned int i = 0; i < len(N1); i++)
        dalphar += N1[i] * i1[i] * std::pow(delta, i1[i]) * std::pow(tau, j1[i]);
    for (unsigned int i = 0; i < len(N2); i++)
        dalphar += N2[i] * std::pow(delta, i2[i]) * std::pow(tau, j2[i]) *
                   std::exp(-std::pow(delta, l2[i])) * (i2[i] - l2[i] * std::pow(delta, l2[i]));
    for (unsigned int i = 0; i < len(N3); i++)
        dalphar += N3[i] * std::pow(delta, i3[i]) * std::pow(tau, j3[i]) *
                   std::exp(-phi3[i] * sqr(delta - 1.0) - beta3[i] * sqr(tau - gamma3[i])) *
                   (i3[i] - 2.0 * delta * phi3[i] * (delta - 1.0));

    return dalpha0 + dalphar / delta;
}

double
Nitrogen::dalpha_dtau(double delta, double tau)
{
    // Ideal gas component of the Helmholtz free energy
    const double dalpha0 = a[0] + a[2] * tau - a[3] / tau - 2.0 * a[4] / sqr(tau) -
                           3.0 * a[5] / cb(tau) + a[6] * a[7] * tau / (std::exp(a[7] * tau) - 1.0);

    // Residual component of the Helmholtz free energy
    double dalphar = 0.0;
    for (unsigned int i = 0; i < len(N1); i++)
        dalphar += N1[i] * j1[i] * std::pow(delta, i1[i]) * std::pow(tau, j1[i]);
    for (unsigned int i = 0; i < len(N2); i++)
        dalphar += N2[i] * j2[i] * std::pow(delta, i2[i]) * std::pow(tau, j2[i]) *
                   std::exp(-std::pow(delta, l2[i]));
    for (unsigned int i = 0; i < len(N3); i++)
        dalphar += N3[i] * std::pow(delta, i3[i]) * std::pow(tau, j3[i]) *
                   std::exp(-phi3[i] * sqr(delta - 1.0) - beta3[i] * sqr(tau - gamma3[i])) *
                   (j3[i] - 2.0 * tau * beta3[i] * (tau - gamma3[i]));

    return (dalpha0 + dalphar) / tau;
}

double
Nitrogen::d2alpha_ddelta2(double delta, double tau)
{
    // Ideal gas component of the Helmholtz free energy
    const double dalpha0 = -1.0 / delta / delta;

    // Residual component of the Helmholtz free energy
    double dalphar = 0.0;
    for (unsigned int i = 0; i < len(N1); i++)
        dalphar += N1[i] * i1[i] * (i1[i] - 1.0) * std::pow(delta, i1[i]) * std::pow(tau, j1[i]);

    for (unsigned int i = 0; i < len(N2); i++)
        dalphar += N2[i] * std::pow(delta, i2[i]) * std::pow(tau, j2[i]) *
                   std::exp(-std::pow(delta, l2[i])) *
                   ((i2[i] - l2[i] * std::pow(delta, l2[i])) *
                        (i2[i] - 1.0 - l2[i] * std::pow(delta, l2[i])) -
                    l2[i] * l2[i] * std::pow(delta, l2[i]));

    for (unsigned int i = 0; i < len(N3); i++)
        dalphar += N3[i] * std::pow(delta, i3[i]) * std::pow(tau, j3[i]) *
                   std::exp(-phi3[i] * sqr(delta - 1.0) - beta3[i] * sqr(tau - gamma3[i])) *
                   (sqr(i3[i] - 2.0 * delta * phi3[i] * (delta - 1.0)) - i3[i] -
                    2.0 * delta * delta * phi3[i]);

    return dalpha0 + dalphar / delta / delta;
}

double
Nitrogen::d2alpha_dtau2(double delta, double tau)
{
    // Ideal gas component of the Helmholtz free energy
    const double dalpha0 =
        -a[0] + 2.0 * a[3] / tau + 6.0 * a[4] / sqr(tau) + 12.0 * a[5] / cb(tau) -
        a[6] * a[7] * a[7] * tau * tau * std::exp(a[7] * tau) / sqr(std::exp(a[7] * tau) - 1.0);

    // Residual component of the Helmholtz free energy
    double dalphar = 0.0;
    for (unsigned int i = 0; i < len(N1); i++)
        dalphar += N1[i] * j1[i] * (j1[i] - 1.0) * std::pow(delta, i1[i]) * std::pow(tau, j1[i]);
    for (unsigned int i = 0; i < len(N2); i++)
        dalphar += N2[i] * j2[i] * (j2[i] - 1.0) * std::pow(delta, i2[i]) * std::pow(tau, j2[i]) *
                   std::exp(-std::pow(delta, l2[i]));
    for (unsigned int i = 0; i < len(N3); i++)
        dalphar += N3[i] * std::pow(delta, i3[i]) * std::pow(tau, j3[i]) *
                   std::exp(-phi3[i] * sqr(delta - 1.0) - beta3[i] * sqr(tau - gamma3[i])) *
                   (sqr(j3[i] - 2.0 * tau * beta3[i] * (tau - gamma3[i])) - j3[i] -
                    2.0 * tau * tau * beta3[i]);

    return (dalpha0 + dalphar) / tau / tau;
}

double
Nitrogen::d2alpha_ddeltatau(double delta, double tau)
{
    // Residual component of the Helmholtz free energy (second derivative of ideal
    // component wrt delta and tau is 0)
    double dalphar = 0.0;
    for (unsigned int i = 0; i < len(N1); i++)
        dalphar += N1[i] * i1[i] * j1[i] * std::pow(delta, i1[i]) * std::pow(tau, j1[i]);
    for (unsigned int i = 0; i < len(N2); i++)
        dalphar += N2[i] * j2[i] * std::pow(delta, i2[i]) * std::pow(tau, j2[i]) *
                   std::exp(-std::pow(delta, l2[i])) * (i2[i] - l2[i] * std::pow(delta, l2[i]));
    for (unsigned int i = 0; i < len(N3); i++)
        dalphar += N3[i] * std::pow(delta, i3[i]) * std::pow(tau, j3[i]) *
                   std::exp(-phi3[i] * sqr(delta - 1.0) - beta3[i] * sqr(tau - gamma3[i])) *
                   (i3[i] - 2.0 * delta * phi3[i] * (delta - 1.0)) *
                   (j3[i] - 2.0 * tau * beta3[i] * (tau - gamma3[i]));

    return dalphar / delta / tau;
}

double
Nitrogen::mu_from_rho_T(double rho, double T)
{
    const double delta = rho / rho_c;
    const double tau = T_c / T;
    const double logTstar = std::log(T / 98.94);

    // The dilute gas component
    double logOmega = 0.0;
    for (unsigned int i = 0; i < len(b_mu); i++)
        logOmega += b_mu[i] * std::pow(logTstar, i);
    const double mu0 =
        0.0266958 * std::sqrt(1000.0 * M * T) / (0.3656 * 0.3656 * std::exp(logOmega));

    // The residual component
    double mur = 0.0;
    for (unsigned int i = 0; i < len(N_mu); i++)
        mur += N_mu[i] * std::pow(tau, t_mu[i]) * std::pow(delta, d_mu[i]) *
               std::exp(-gamma_mu[i] * std::pow(delta, l_mu[i]));

    // [Pa.s]
    return (mu0 + mur) * 1.0e-6;
}

double
Nitrogen::k_from_rho_T(double rho, double T)
{
    const double delta = rho / rho_c;
    const double tau = T_c / T;

    // The dilute gas component
    const double lambda0 =
        1.511 * mu_from_rho_T(0.0, T) * 1.0e6 + 2.117 / tau - 3.332 * std::pow(tau, -0.7);

    // The residual component
    double lambdar = 0.0;
    for (unsigned int i = 0; i < len(N_k); ++i)
        lambdar += N_k[i] * std::pow(tau, t_k[i]) * std::pow(delta, d_k[i]) *
                   std::exp(-gamma_k[i] * std::pow(delta, l_k[i]));

    // The thermal conductivity (note: critical enhancement not implemented)
    return (lambda0 + lambdar) * 1.0e-3;
}

} // namespace fprops
