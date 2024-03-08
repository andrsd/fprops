#include "fprops/IdealGas.h"
#include "fprops/Exception.h"
#include <cmath>

namespace fprops {

/// Universal gas constant
static const double R = 8.3144598;

IdealGas::IdealGas(double gamma, double molar_mass) :
    SinglePhaseFluidProperties(),
    gamma(gamma),
    molar_mass(molar_mass),
    R_specific(R / molar_mass),
    cp(gamma * R_specific / (gamma - 1.0)),
    cv(cp / gamma),
    mu(0.),
    k(0.)
{
}

double
IdealGas::get_gamma() const
{
    return this->gamma;
}

double
IdealGas::get_specific_gas_constant() const
{
    return this->R_specific;
}

State
IdealGas::rho_T(double rho, double T) const
{
    if (rho < 0)
        throw Exception("Negative density");
    if (T < 0)
        throw Exception("Negative temperature");

    State state;
    state.rho = rho;
    state.T = T;
    state.cp = this->cp;
    state.cv = this->cv;
    state.mu = this->mu;
    state.k = this->k;
    state.p = rho * R * T / this->molar_mass;
    state.u = this->cv * T;
    state.v = 1. / state.rho;
    const double n = std::pow(T, this->gamma) / std::pow(state.p, this->gamma - 1.0);
    if (n <= 0)
        throw Exception("Invalid log base for computing entropy");
    state.s = this->cv * std::log(n);
    state.h = this->cp * T;
    state.w = std::sqrt(this->cp * R * T / (this->cv * this->molar_mass));
    return state;
}

State
IdealGas::rho_p(double rho, double p) const
{
    if (rho < 0)
        throw Exception("Negative density");

    State state;
    state.rho = rho;
    state.p = p;
    state.cp = this->cp;
    state.cv = this->cv;
    state.mu = this->mu;
    state.k = this->k;
    state.T = p * this->molar_mass / (rho * R);
    state.u = this->cv * state.T;
    state.v = 1. / state.rho;
    const double n = std::pow(state.T, this->gamma) / std::pow(state.p, this->gamma - 1.0);
    if (n <= 0)
        throw Exception("Invalid log base for computing entropy");
    state.s = this->cv * std::log(n);
    state.h = this->cp * state.T;
    state.w = std::sqrt(this->cp * R * state.T / (this->cv * this->molar_mass));
    return state;
}

State
IdealGas::p_T(double p, double T) const
{
    if (T < 0)
        throw Exception("Negative temperature");

    State state;
    state.p = p;
    state.T = T;
    state.cp = this->cp;
    state.cv = this->cv;
    state.mu = this->mu;
    state.k = this->k;
    state.rho = p * this->molar_mass / (R * T);
    state.u = this->cv * T;
    state.v = 1. / state.rho;
    const double n = std::pow(T, this->gamma) / std::pow(p, this->gamma - 1.0);
    if (n <= 0)
        throw Exception("Invalid log base for computing entropy");
    state.s = this->cv * std::log(n);
    state.h = this->cp * T;
    state.w = std::sqrt(this->cp * R * T / (this->cv * this->molar_mass));
    return state;
}

State
IdealGas::v_u(double v, double u) const
{
    if (v <= 0.)
        throw Exception("Negative specific volume");
    if (u <= 0.)
        throw Exception("Negative internal energy");

    State state;
    state.v = v;
    state.u = u;
    state.cp = this->cp;
    state.cv = this->cv;
    state.mu = this->mu;
    state.k = this->k;
    state.rho = 1. / v;
    state.p = (this->gamma - 1.0) * u * state.rho;
    state.T = u / this->cv;
    const double n = std::pow(state.T, this->gamma) / std::pow(state.p, this->gamma - 1.0);
    if (n <= 0)
        throw Exception("Invalid log base for computing entropy");
    state.s = this->cv * std::log(n);
    state.h = this->cp * state.T;
    state.w = std::sqrt(this->gamma * this->R_specific * state.T);
    return state;
}

State
IdealGas::h_s(double h, double s) const
{
    State state;
    state.h = h;
    state.s = s;
    state.cp = this->cp;
    state.cv = this->cv;
    state.mu = this->mu;
    state.k = this->k;
    state.p = std::pow(h / (this->gamma * this->cv), this->gamma / (this->gamma - 1.0)) *
              std::exp(-s / ((this->gamma - 1.0) * this->cv));
    const double aux = (s + this->cv * std::log(std::pow(state.p, this->gamma - 1.0))) / this->cv;
    state.T = std::pow(std::exp(aux), 1.0 / this->gamma);
    state.rho = state.p * this->molar_mass / (R * state.T);
    state.u = this->cv * state.T;
    state.v = 1. / state.rho;
    const double n = std::pow(state.T, this->gamma) / std::pow(state.p, this->gamma - 1.0);
    if (n <= 0)
        throw Exception("Invalid log base for computing entropy");
    state.w = std::sqrt(this->gamma * this->R_specific * state.T);
    return state;
}

void
IdealGas::set_mu(double mu)
{
    this->mu = mu;
}

void
IdealGas::set_k(double k)
{
    this->k = k;
}

} // namespace fprops
