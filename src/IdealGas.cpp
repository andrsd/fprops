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
