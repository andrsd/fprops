#include "IdealGas.h"
#include <cmath>
#include <stdexcept>

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

SinglePhaseFluidProperties::Props
IdealGas::rho_T(double rho, double T) const
{
    if (rho < 0)
        throw std::domain_error("Negative density");
    if (T < 0)
        throw std::domain_error("Negative temperature");

    Props props;
    props.rho = rho;
    props.T = T;
    props.cp = this->cp;
    props.cv = this->cv;
    props.mu = this->mu;
    props.k = this->k;
    props.p = rho * R * T / this->molar_mass;
    props.u = this->cv * T;
    props.v = 1. / props.rho;
    const double n = std::pow(T, this->gamma) / std::pow(props.p, this->gamma - 1.0);
    if (n <= 0)
        throw std::domain_error("Invalid log base for computing entropy");
    props.s = this->cv * std::log(n);
    props.h = this->cp * T;
    props.w = std::sqrt(this->cp * R * T / (this->cv * this->molar_mass));
    return props;
}

SinglePhaseFluidProperties::Props
IdealGas::rho_p(double rho, double p) const
{
    if (rho < 0)
        throw std::domain_error("Negative density");

    Props props;
    props.rho = rho;
    props.p = p;
    props.cp = this->cp;
    props.cv = this->cv;
    props.mu = this->mu;
    props.k = this->k;
    props.T = p * this->molar_mass / (rho * R);
    props.u = this->cv * props.T;
    props.v = 1. / props.rho;
    const double n = std::pow(props.T, this->gamma) / std::pow(props.p, this->gamma - 1.0);
    if (n <= 0)
        throw std::domain_error("Invalid log base for computing entropy");
    props.s = this->cv * std::log(n);
    props.h = this->cp * props.T;
    props.w = std::sqrt(this->cp * R * props.T / (this->cv * this->molar_mass));
    return props;
}

SinglePhaseFluidProperties::Props
IdealGas::p_T(double p, double T) const
{
    if (T < 0)
        throw std::domain_error("Negative temperature");

    Props props;
    props.p = p;
    props.T = T;
    props.cp = this->cp;
    props.cv = this->cv;
    props.mu = this->mu;
    props.k = this->k;
    props.rho = p * this->molar_mass / (R * T);
    props.u = this->cv * T;
    props.v = 1. / props.rho;
    const double n = std::pow(T, this->gamma) / std::pow(p, this->gamma - 1.0);
    if (n <= 0)
        throw std::domain_error("Invalid log base for computing entropy");
    props.s = this->cv * std::log(n);
    props.h = this->cp * T;
    props.w = std::sqrt(this->cp * R * T / (this->cv * this->molar_mass));
    return props;
}

SinglePhaseFluidProperties::Props
IdealGas::v_u(double v, double u) const
{
    if (v <= 0.)
        throw std::domain_error("Negative specific volume");
    if (u <= 0.)
        throw std::domain_error("Negative internal energy");

    Props props;
    props.v = v;
    props.u = u;
    props.cp = this->cp;
    props.cv = this->cv;
    props.mu = this->mu;
    props.k = this->k;
    props.rho = 1. / v;
    props.p = (this->gamma - 1.0) * u * props.rho;
    props.T = u / this->cv;
    const double n = std::pow(props.T, this->gamma) / std::pow(props.p, this->gamma - 1.0);
    if (n <= 0)
        throw std::domain_error("Invalid log base for computing entropy");
    props.s = this->cv * std::log(n);
    props.h = this->cp * props.T;
    props.w = std::sqrt(this->gamma * this->R_specific * props.T);
    return props;
}

SinglePhaseFluidProperties::Props
IdealGas::h_s(double h, double s) const
{
    Props props;
    props.h = h;
    props.s = s;
    props.cp = this->cp;
    props.cv = this->cv;
    props.mu = this->mu;
    props.k = this->k;
    props.p = std::pow(h / (this->gamma * this->cv), this->gamma / (this->gamma - 1.0)) *
              std::exp(-s / ((this->gamma - 1.0) * this->cv));
    const double aux = (s + this->cv * std::log(std::pow(props.p, this->gamma - 1.0))) / this->cv;
    props.T = std::pow(std::exp(aux), 1.0 / this->gamma);
    props.rho = props.p * this->molar_mass / (R * props.T);
    props.u = this->cv * props.T;
    props.v = 1. / props.rho;
    const double n = std::pow(props.T, this->gamma) / std::pow(props.p, this->gamma - 1.0);
    if (n <= 0)
        throw std::domain_error("Invalid log base for computing entropy");
    props.w = std::sqrt(this->gamma * this->R_specific * props.T);
    return props;
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
