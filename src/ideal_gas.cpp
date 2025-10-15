// SPDX-FileCopyrightText: 2022 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#include "fprops/ideal_gas.h"
#include "fprops/exception.h"
#include <cmath>

namespace fprops {

/// Universal gas constant
static const double R = 8.3144598;

IdealGas::IdealGas(double gamma, double molar_mass) :
    m_gamma(gamma),
    m_molar_mass(molar_mass),
    m_R_specific(R / molar_mass),
    m_cp(gamma * m_R_specific / (gamma - 1.0)),
    m_cv(m_cp / gamma),
    m_mu(0.),
    m_k(0.)
{
}

double
IdealGas::gamma() const
{
    return this->m_gamma;
}

double
IdealGas::molar_mass() const
{
    return this->m_molar_mass;
}

double
IdealGas::R_specific() const
{
    return this->m_R_specific;
}

double
IdealGas::cp() const
{
    return this->m_cp;
}

double
IdealGas::cv() const
{
    return this->m_cv;
}

double
IdealGas::mu() const
{
    return this->m_mu;
}

double
IdealGas::k() const
{
    return this->m_k;
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
    state.cp = this->m_cp;
    state.cv = this->m_cv;
    state.mu = this->m_mu;
    state.k = this->m_k;
    state.p = rho * R * T / this->m_molar_mass;
    state.u = this->m_cv * T;
    state.v = 1. / state.rho;
    const double n = std::pow(T, this->m_gamma) / std::pow(state.p, this->m_gamma - 1.0);
    if (n <= 0)
        throw Exception("Invalid log base for computing entropy");
    state.s = this->m_cv * std::log(n);
    state.h = this->m_cp * T;
    state.w = std::sqrt(this->m_cp * R * T / (this->m_cv * this->m_molar_mass));
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
    state.cp = this->m_cp;
    state.cv = this->m_cv;
    state.mu = this->m_mu;
    state.k = this->m_k;
    state.T = p * this->m_molar_mass / (rho * R);
    state.u = this->m_cv * state.T;
    state.v = 1. / state.rho;
    const double n = std::pow(state.T, this->m_gamma) / std::pow(state.p, this->m_gamma - 1.0);
    if (n <= 0)
        throw Exception("Invalid log base for computing entropy");
    state.s = this->m_cv * std::log(n);
    state.h = this->m_cp * state.T;
    state.w = std::sqrt(this->m_cp * R * state.T / (this->m_cv * this->m_molar_mass));
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
    state.cp = this->m_cp;
    state.cv = this->m_cv;
    state.mu = this->m_mu;
    state.k = this->m_k;
    state.rho = p * this->m_molar_mass / (R * T);
    state.u = this->m_cv * T;
    state.v = 1. / state.rho;
    const double n = std::pow(T, this->m_gamma) / std::pow(p, this->m_gamma - 1.0);
    if (n <= 0)
        throw Exception("Invalid log base for computing entropy");
    state.s = this->m_cv * std::log(n);
    state.h = this->m_cp * T;
    state.w = std::sqrt(this->m_cp * R * T / (this->m_cv * this->m_molar_mass));
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
    state.cp = this->m_cp;
    state.cv = this->m_cv;
    state.mu = this->m_mu;
    state.k = this->m_k;
    state.rho = 1. / v;
    state.p = (this->m_gamma - 1.0) * u * state.rho;
    state.T = u / this->m_cv;
    const double n = std::pow(state.T, this->m_gamma) / std::pow(state.p, this->m_gamma - 1.0);
    if (n <= 0)
        throw Exception("Invalid log base for computing entropy");
    state.s = this->m_cv * std::log(n);
    state.h = this->m_cp * state.T;
    state.w = std::sqrt(this->m_gamma * this->m_R_specific * state.T);
    return state;
}

State
IdealGas::h_s(double h, double s) const
{
    State state;
    state.h = h;
    state.s = s;
    state.cp = this->m_cp;
    state.cv = this->m_cv;
    state.mu = this->m_mu;
    state.k = this->m_k;
    state.p = std::pow(h / (this->m_gamma * this->m_cv), this->m_gamma / (this->m_gamma - 1.0)) *
              std::exp(-s / ((this->m_gamma - 1.0) * this->m_cv));
    const double aux =
        (s + this->m_cv * std::log(std::pow(state.p, this->m_gamma - 1.0))) / this->m_cv;
    state.T = std::pow(std::exp(aux), 1.0 / this->m_gamma);
    state.rho = state.p * this->m_molar_mass / (R * state.T);
    state.u = this->m_cv * state.T;
    state.v = 1. / state.rho;
    const double n = std::pow(state.T, this->m_gamma) / std::pow(state.p, this->m_gamma - 1.0);
    if (n <= 0)
        throw Exception("Invalid log base for computing entropy");
    state.w = std::sqrt(this->m_gamma * this->m_R_specific * state.T);
    return state;
}

State
IdealGas::v_h(double v, double h) const
{
    if (v < 0)
        throw Exception("Negative specific volume");

    State state;
    state.v = v;
    state.h = h;
    state.cp = this->m_cp;
    state.cv = this->m_cv;
    state.mu = this->m_mu;
    state.k = this->m_k;
    state.rho = 1. / v;
    state.T = h / this->m_cp;
    state.p = state.rho * R * state.T / this->m_molar_mass;
    state.u = this->m_cv * state.T;
    const double n = std::pow(state.T, this->m_gamma) / std::pow(state.p, this->m_gamma - 1.0);
    if (n <= 0)
        throw Exception("Invalid log base for computing entropy");
    state.s = this->m_cv * std::log(n);
    state.w = std::sqrt(this->m_cp * R * state.T / (this->m_cv * this->m_molar_mass));
    return state;
}

void
IdealGas::set_mu(double mu)
{
    this->m_mu = mu;
}

void
IdealGas::set_k(double k)
{
    this->m_k = k;
}

} // namespace fprops
