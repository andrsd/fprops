#include "Helmholtz.h"
#include "BrentsMethod.h"
#include "Numerics.h"
#include <cmath>
#include <assert.h>

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
Helmholtz::p_T(double p, double T)
{
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

    // u
    const double u = this->R * T * tau * da_dt / this->M;

    // h
    const double h = this->R * T * (tau * da_dt + delta * da_dd) / this->M;

    // w
    const double n = 2.0 * delta * da_dd + delta * delta * d2a_dd2 -
                     sqr(delta * da_dd - delta * tau * d2a_ddt) / (tau * tau * d2a_dt2);
    const double w = std::sqrt(this->R * T * n / this->M);

    // cp = dh/dt
    const double cp = this->R *
                      (-tau * tau * d2a_dt2 + sqr(delta * da_dd - delta * tau * d2a_ddt) /
                                                  (2.0 * delta * da_dd + delta * delta * d2a_dd2)) /
                      this->M;

    // cv = du/dt
    const double cv = -this->R * tau * tau * d2a_dt2 / this->M;

    // s = ...
    const double s = this->R * (tau * da_dt - a) / this->M;

    // mu
    const double mu = mu_from_rho_T(rho, T);

    // k
    const double k = k_from_rho_T(rho, T);

    props.p = p;
    props.T = T;
    props.cp = cp;
    props.cv = cv;
    props.mu = mu;
    props.k = k;
    props.rho = rho;
    props.u = u;
    props.v = 1. / rho;
    props.s = s;
    props.h = h;
    props.w = w;

    return props;
}

SinglePhaseFluidProperties::Props
Helmholtz::v_u(double v, double u)
{
    assert(v != 0.);

    Props props;
    props.v = v;
    props.u = u;
    return props;
}

double
Helmholtz::p_from_rho_T(double density, double temperature)
{
    // Scale the input density and temperature
    const double delta = density / this->rho_c;
    const double tau = this->T_c / temperature;

    return this->R * density * temperature * delta * dalpha_ddelta(delta, tau) / this->M;
}

double
Helmholtz::rho_from_p_T(double pressure, double temperature)
{
    double density;
    // Initial estimate of a bracketing interval for the density
    double lower_density = 1.0e-2;
    double upper_density = 100.0;

    // The density is found by finding the zero of the pressure
    auto pressure_diff = [&pressure, &temperature, this](double x) {
        return this->p_from_rho_T(x, temperature) - pressure;
    };

    BrentsMethod::bracket(pressure_diff, lower_density, upper_density);
    density = BrentsMethod::root(pressure_diff, lower_density, upper_density);

    return density;
}

} // namespace fprops
