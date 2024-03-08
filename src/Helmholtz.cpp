#include "fprops/Helmholtz.h"
#include "fprops/Numerics.h"
#include "fprops/Exception.h"
#include <cmath>

namespace fprops {

Helmholtz::Helmholtz(double R, double M, double rho_c, double T_c) :
    SinglePhaseFluidProperties(),
    R(R),
    M(M),
    rho_c(rho_c),
    T_c(T_c)
{
}

State
Helmholtz::rho_T(double rho, double T) const
{
    if (rho < 0)
        throw Exception("Negative density");
    if (T < 0)
        throw Exception("Negative temperature");

    const double delta = rho / this->rho_c;
    const double tau = this->T_c / T;

    const double a = alpha(delta, tau);
    const double da_dd = dalpha_ddelta(delta, tau);
    const double da_dt = dalpha_dtau(delta, tau);
    const double d2a_dt2 = d2alpha_dtau2(delta, tau);
    const double d2a_dd2 = d2alpha_ddelta2(delta, tau);
    const double d2a_ddt = d2alpha_ddeltatau(delta, tau);

    auto v = 1. / rho;
    auto p = pressure(rho, T, delta, da_dd);
    auto u = internal_energy(T, tau, da_dt);
    auto h = enthalphy(T, delta, tau, da_dt, da_dd);
    auto w = sound_speed(T, delta, tau, da_dd, d2a_dd2, d2a_ddt, d2a_dt2);
    auto cp = heat_capacity_isobaric(delta, tau, da_dd, d2a_dt2, d2a_dd2, d2a_ddt);
    auto cv = heat_capacity_isochoric(tau, d2a_dt2);
    auto s = entropy(tau, a, da_dt);
    auto mu = mu_from_rho_T(rho, T);
    auto k = k_from_rho_T(rho, T);

    return State(u, v, rho, p, T, mu, cp, cv, s, k, h, w);
}

State
Helmholtz::rho_p(double rho, double p) const
{
    if (rho < 0)
        throw Exception("Negative density");

    const double T = T_from_rho_p(rho, p);

    const double delta = rho / this->rho_c;
    const double tau = this->T_c / T;

    const double a = alpha(delta, tau);
    const double da_dd = dalpha_ddelta(delta, tau);
    const double da_dt = dalpha_dtau(delta, tau);
    const double d2a_dt2 = d2alpha_dtau2(delta, tau);
    const double d2a_dd2 = d2alpha_ddelta2(delta, tau);
    const double d2a_ddt = d2alpha_ddeltatau(delta, tau);

    auto v = 1. / rho;
    auto u = internal_energy(T, tau, da_dt);
    auto h = enthalphy(T, delta, tau, da_dt, da_dd);
    auto w = sound_speed(T, delta, tau, da_dd, d2a_dd2, d2a_ddt, d2a_dt2);
    auto cp = heat_capacity_isobaric(delta, tau, da_dd, d2a_dt2, d2a_dd2, d2a_ddt);
    auto cv = heat_capacity_isochoric(tau, d2a_dt2);
    auto s = entropy(tau, a, da_dt);
    auto mu = mu_from_rho_T(rho, T);
    auto k = k_from_rho_T(rho, T);

    return State(u, v, rho, p, T, mu, cp, cv, s, k, h, w);
}

State
Helmholtz::p_T(double p, double T) const
{
    if (T < 0)
        throw Exception("Negative temperature");

    const double rho = rho_from_p_T(p, T);

    const double delta = rho / this->rho_c;
    const double tau = this->T_c / T;

    const double a = alpha(delta, tau);
    const double da_dd = dalpha_ddelta(delta, tau);
    const double da_dt = dalpha_dtau(delta, tau);
    const double d2a_dt2 = d2alpha_dtau2(delta, tau);
    const double d2a_dd2 = d2alpha_ddelta2(delta, tau);
    const double d2a_ddt = d2alpha_ddeltatau(delta, tau);

    auto v = 1. / rho;
    auto u = internal_energy(T, tau, da_dt);
    auto h = enthalphy(T, delta, tau, da_dt, da_dd);
    auto w = sound_speed(T, delta, tau, da_dd, d2a_dd2, d2a_ddt, d2a_dt2);
    auto cp = heat_capacity_isobaric(delta, tau, da_dd, d2a_dt2, d2a_dd2, d2a_ddt);
    auto cv = heat_capacity_isochoric(tau, d2a_dt2);
    auto s = entropy(tau, a, da_dt);
    auto mu = mu_from_rho_T(rho, T);
    auto k = k_from_rho_T(rho, T);

    return State(u, v, rho, p, T, mu, cp, cv, s, k, h, w);
}

State
Helmholtz::v_u(double v, double u) const
{
    if (v <= 0.)
        throw Exception("Negative specific volume");
    if (u <= 0.)
        throw Exception("Negative internal energy");

    const double rho = 1. / v;
    const double delta = rho / this->rho_c;
    const double tau = tau_from_v_u(v, u);

    const double a = alpha(delta, tau);
    const double da_dd = dalpha_ddelta(delta, tau);
    const double da_dt = dalpha_dtau(delta, tau);
    const double d2a_dt2 = d2alpha_dtau2(delta, tau);
    const double d2a_dd2 = d2alpha_ddelta2(delta, tau);
    const double d2a_ddt = d2alpha_ddeltatau(delta, tau);

    auto T = temperature(u, tau, da_dt);
    auto p = pressure(rho, T, delta, da_dd);
    auto h = enthalphy(T, delta, tau, da_dt, da_dd);
    auto w = sound_speed(T, delta, tau, da_dd, d2a_dd2, d2a_ddt, d2a_dt2);
    auto cp = heat_capacity_isobaric(delta, tau, da_dd, d2a_dt2, d2a_dd2, d2a_ddt);
    auto cv = heat_capacity_isochoric(tau, d2a_dt2);
    auto s = entropy(tau, a, da_dt);
    auto mu = mu_from_rho_T(rho, T);
    auto k = k_from_rho_T(rho, T);

    return State(u, v, rho, p, T, mu, cp, cv, s, k, h, w);
}

State
Helmholtz::h_s(double h, double s) const
{
    throw Exception("Not implemented");
}

double
Helmholtz::rho_from_p_T(double p, double T) const
{
    auto f = [&p, &T, this](double rho) {
        const double delta = rho / this->rho_c;
        const double tau = this->T_c / T;

        return this->R * rho * T * delta * dalpha_ddelta(delta, tau) / this->M - p;
    };
    auto df = [&T, this](double rho) {
        const double delta = rho / this->rho_c;
        const double ddelta_drho = 1 / this->rho_c;
        const double tau = this->T_c / T;

        double K = this->R * T / this->M;
        double t1 = delta * dalpha_ddelta(delta, tau);
        double t2 = rho * ddelta_drho * dalpha_ddelta(delta, tau);
        double t3 = rho * delta * d2alpha_ddelta2(delta, tau) * ddelta_drho;

        return K * (t1 + t2 + t3);
    };

    return newton::root(1.0e-2, f, df);
}

double
Helmholtz::T_from_rho_p(double rho, double p) const
{
    auto f = [&rho, &p, this](double T) {
        const double delta = rho / this->rho_c;
        const double tau = this->T_c / T;

        return this->R * rho * T * delta * dalpha_ddelta(delta, tau) / this->M - p;
    };
    auto df = [&rho, this](double T) {
        const double delta = rho / this->rho_c;
        const double tau = this->T_c / T;
        const double dtau_dT = -this->T_c / T / T;

        double K = this->R * rho * delta / this->M;
        double t1 = dalpha_ddelta(delta, tau);
        double t2 = T * d2alpha_ddeltatau(delta, tau) * dtau_dT;

        return K * (t1 + t2);
    };

    return newton::root(275., f, df);
}

double
Helmholtz::tau_from_v_u(double v, double u) const
{
    auto f = [&v, &u, this](double tau) {
        const double delta = 1. / this->rho_c / v;
        return this->R * this->T_c * dalpha_dtau(delta, tau) / this->M - u;
    };
    auto df = [&v, this](double tau) {
        const double delta = 1. / this->rho_c / v;
        return this->R * this->T_c * d2alpha_dtau2(delta, tau) / this->M;
    };

    return newton::root(1e-1, f, df);
}

double
Helmholtz::temperature(double u, double tau, double da_dt) const
{
    return u * this->M / (this->R * tau * da_dt);
}

double
Helmholtz::pressure(double rho, double T, double delta, double da_dd) const
{
    return this->R * rho * T * delta * da_dd / this->M;
}

double
Helmholtz::internal_energy(double T, double tau, double da_dt) const
{
    return this->R * T * tau * da_dt / this->M;
}

double
Helmholtz::enthalphy(double T, double delta, double tau, double da_dt, double da_dd) const
{
    return this->R * T * (tau * da_dt + delta * da_dd) / this->M;
}

double
Helmholtz::sound_speed(double T,
                       double delta,
                       double tau,
                       double da_dd,
                       double d2a_dd2,
                       double d2a_ddt,
                       double d2a_dt2) const
{
    const double n = 2.0 * delta * da_dd + delta * delta * d2a_dd2 -
                     sqr(delta * da_dd - delta * tau * d2a_ddt) / (tau * tau * d2a_dt2);
    return std::sqrt(this->R * T * n / this->M);
}

double
Helmholtz::entropy(double tau, double a, double da_dt) const
{
    return this->R * (tau * da_dt - a) / this->M;
}

double
Helmholtz::heat_capacity_isobaric(double delta,
                                  double tau,
                                  double da_dd,
                                  double d2a_dt2,
                                  double d2a_dd2,
                                  double d2a_ddt) const
{
    return this->R *
           (-tau * tau * d2a_dt2 + sqr(delta * da_dd - delta * tau * d2a_ddt) /
                                       (2.0 * delta * da_dd + delta * delta * d2a_dd2)) /
           this->M;
}

double
Helmholtz::heat_capacity_isochoric(double tau, double d2a_dt2) const
{
    return -this->R * tau * tau * d2a_dt2 / this->M;
}

} // namespace fprops
