#include "Helmholtz.h"
#include "Numerics.h"
#include <cmath>
#include <stdexcept>

namespace fprops {

Helmholtz::Helmholtz(double R, double M, double rho_c, double T_c) :
    SinglePhaseFluidProperties(),
    R(R),
    M(M),
    rho_c(rho_c),
    T_c(T_c)
{
}

SinglePhaseFluidProperties::Props
Helmholtz::p_T(double p, double T) const
{
    if (T < 0)
        throw std::domain_error("Negative temperature");

    Props props;

    const double rho = rho_from_p_T(p, T);

    const double delta = rho / this->rho_c;
    const double tau = this->T_c / T;

    const double a = alpha(delta, tau);
    const double da_dd = dalpha_ddelta(delta, tau);
    const double da_dt = dalpha_dtau(delta, tau);
    const double d2a_dt2 = d2alpha_dtau2(delta, tau);
    const double d2a_dd2 = d2alpha_ddelta2(delta, tau);
    const double d2a_ddt = d2alpha_ddeltatau(delta, tau);

    props.p = p;
    props.T = T;
    props.rho = rho;
    props.v = 1. / rho;
    // u
    props.u = this->R * T * tau * da_dt / this->M;
    // h
    props.h = this->R * T * (tau * da_dt + delta * da_dd) / this->M;
    // w
    const double n = 2.0 * delta * da_dd + delta * delta * d2a_dd2 -
                     sqr(delta * da_dd - delta * tau * d2a_ddt) / (tau * tau * d2a_dt2);
    props.w = std::sqrt(this->R * T * n / this->M);
    // cp = dh/dt
    props.cp = this->R *
               (-tau * tau * d2a_dt2 + sqr(delta * da_dd - delta * tau * d2a_ddt) /
                                           (2.0 * delta * da_dd + delta * delta * d2a_dd2)) /
               this->M;
    // cv = du/dt
    props.cv = -this->R * tau * tau * d2a_dt2 / this->M;
    // s = ...
    props.s = this->R * (tau * da_dt - a) / this->M;
    // mu
    props.mu = mu_from_rho_T(rho, T);
    // k
    props.k = k_from_rho_T(rho, T);
    return props;
}

SinglePhaseFluidProperties::Props
Helmholtz::v_u(double v, double u) const
{
    if (v <= 0.)
        throw std::domain_error("Negative specific volume");
    if (u <= 0.)
        throw std::domain_error("Negative internal energy");

    Props props;

    const double rho = 1. / v;
    const double delta = rho / this->rho_c;
    const double tau = tau_from_v_u(v, u);

    const double a = alpha(delta, tau);
    const double da_dd = dalpha_ddelta(delta, tau);
    const double da_dt = dalpha_dtau(delta, tau);
    const double d2a_dt2 = d2alpha_dtau2(delta, tau);
    const double d2a_dd2 = d2alpha_ddelta2(delta, tau);
    const double d2a_ddt = d2alpha_ddeltatau(delta, tau);

    props.rho = rho;
    props.u = u;
    props.v = v;
    // T
    props.T = u * this->M / (this->R * tau * da_dt);
    // p
    props.p = this->R * rho * props.T * delta * da_dd / this->M;
    // h
    props.h = this->R * props.T * (tau * da_dt + delta * da_dd) / this->M;
    // w
    const double n = 2.0 * delta * da_dd + delta * delta * d2a_dd2 -
                     sqr(delta * da_dd - delta * tau * d2a_ddt) / (tau * tau * d2a_dt2);
    props.w = std::sqrt(this->R * props.T * n / this->M);
    // cp = dh/dt
    props.cp = this->R *
               (-tau * tau * d2a_dt2 + sqr(delta * da_dd - delta * tau * d2a_ddt) /
                                           (2.0 * delta * da_dd + delta * delta * d2a_dd2)) /
               this->M;
    // cv = du/dt
    props.cv = -this->R * tau * tau * d2a_dt2 / this->M;
    // s = ...
    props.s = this->R * (tau * da_dt - a) / this->M;
    // mu
    props.mu = mu_from_rho_T(rho, props.T);
    // k
    props.k = k_from_rho_T(rho, props.T);
    return props;
}

SinglePhaseFluidProperties::Props
Helmholtz::h_s(double h, double s) const
{
    throw std::domain_error("Not implemented");
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

} // namespace fprops
