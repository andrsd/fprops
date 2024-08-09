// SPDX-FileCopyrightText: 2022 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#pragma once

#include "fprops/numerics.h"
#include <vector>
#include <cassert>

namespace fprops {

/// Model for computing viscosity (eta_0)
///
/// @tparam T The basic data type
///
/// \f$ lambda_0 = A_0 * \eta_0 + \displaystyle\sum_{i=1}^{n} A_i \cdot \tau^{t_i} \f$
template <typename T>
class Eta0AndPoly {
public:
    /// Eta0 and polynomial model
    ///
    /// @param A \f$A_i\f$ coefficients
    /// @param t \f$t_i\f$ coefficients
    ///
    /// Note that coefficient \f$t_0\f$ is not used
    Eta0AndPoly(const std::vector<double> & A, const std::vector<double> & t) : A(A), t(t)
    {
        assert(A.size() >= 1);
        assert(A.size() == t.size());
    }

    /// Evaluate the model
    ///
    /// @param eta0 \f$\eta_0\f$
    /// @param tau  \f$\tau\f$
    /// @return The computed value
    T
    value(double eta0, double tau) const
    {
        T sum = this->A[0] * eta0;
        for (unsigned int i = 1; i < A.size(); i++)
            sum += this->A[i] * math::pow(tau, this->t[i]);
        return sum;
    }

private:
    std::vector<double> A;
    std::vector<double> t;
};

/// Lennard-Jones model for computing viscosity
///
/// @tparam T The basic data type
///
/// \f$ \eta_0(T) = \frac{C \sqrt{M T}}{\sigma^2 \Omega(T^*)} \f$
///
/// where \f$\sigma\f$ is the Lennard-Jones size parameter and
/// \f$\Omega\f$ is the collision integral, given by
///
/// \f$ \Omega(T^*) = \exp(\displaystyle\sum_{i=0}^n b_i \ln(T^*))^i \f$
///
/// where \f$T^* = T / (\epsilon / k)) \f$ and \f$\epsilon / k\f$ is the Lennard-Jones
/// energy parameter.
template <typename T>
class LennardJones {
public:
    /// Lennard-Jones model
    ///
    /// @param C Constant in front of term
    /// @param M Molar mass \f$[{kg\over mol}]\f$
    /// @param epsilon_over_k \f${\epsilon\over k} [K]\f$
    /// @param sigma \f$\sigma\f$
    /// @param b \f$b_i\f$
    LennardJones(double C,
                 double M,
                 double epsilon_over_k,
                 double sigma,
                 const std::vector<double> & b) :
        C(C),
        M(M),
        epsilon_over_k(epsilon_over_k),
        sigma(sigma),
        b(b)
    {
    }

    /// Evaluate the model
    ///
    /// @param temperature Temperature \f$[K]\f$
    /// @return The computed value
    T
    value(double temperature) const
    {
        double log_T_star = std::log(temperature / this->epsilon_over_k);
        double Omega_T_star = 0;
        for (unsigned int i = 0; i < this->b.size(); i++)
            Omega_T_star += this->b[i] * math::pow(log_T_star, i);
        Omega_T_star = std::exp(Omega_T_star);
        return this->C * std::sqrt(1000.0 * this->M * temperature) /
               (this->sigma * this->sigma * Omega_T_star);
    }

private:
    double C;
    double M;
    double epsilon_over_k;
    const double sigma;
    std::vector<double> b;
};

/// Modified Batshinski-Hildebrand model
///
/// @tparam T The basic data type
///
/// \f$ v = \displaystyle\sum_{i=0}^{n} N_i \tau^{t_i} \delta^{d_i} \exp(-\gamma_i \delta^{l_i})\f$
template <typename T>
class ModifiedBatshinskiHildebrand {
public:
    /// Modified Batshinki-Hildebrand
    ///
    /// @param n \f$N_i\f$
    /// @param t \f$t_i\f$
    /// @param d \f$d_i\f$
    /// @param gamma \f$\gamma_i\f$
    /// @param l \f$l_i\f$
    ModifiedBatshinskiHildebrand(const std::vector<double> & n,
                                 const std::vector<double> & t,
                                 const std::vector<double> & d,
                                 const std::vector<double> & gamma,
                                 const std::vector<double> & l) :
        n(n),
        t(t),
        d(d),
        gamma(gamma),
        l(l)
    {
    }

    /// Evaluate the model
    ///
    /// @param delta \f$\delta\f$
    /// @param tau \f$\tau\f$
    /// @return Computed value
    T
    value(double delta, double tau) const
    {
        double sum = 0.0;
        for (unsigned int i = 0; i < n.size(); i++)
            sum += this->n[i] * math::pow(tau, this->t[i]) * math::pow(delta, this->d[i]) *
                   std::exp(-this->gamma[i] * math::pow(delta, this->l[i]));
        return sum;
    }

private:
    std::vector<double> n;
    std::vector<double> t;
    std::vector<double> d;
    std::vector<double> gamma;
    std::vector<double> l;
};

/// Powers of temperature
///
/// @tparam T The basic data type
///
/// \f$ v = \displaystyle\sum_{i=0}^{n} a_i T^{t_i}\f$
template <typename TYPE>
class PowersOfTemperature {
public:
    ///
    ///
    /// @param a Array of \f$a_i\f$ coefficients
    /// @param t Array of \f$t_i\f$ exponents
    PowersOfTemperature(const std::vector<double> & a, const std::vector<double> & t) : a(a), t(t)
    {
    }

    /// Evaluate the model
    ///
    /// @param T Temperature
    /// @return Computed value
    TYPE
    value(TYPE T) const
    {
        TYPE sum = 0;
        for (std::size_t i = 0; i < a.size(); ++i)
            sum += this->a[i] * math::pow(T, this->t[i]);
        return sum;
    }

private:
    /// a_i coefficients
    std::vector<double> a;
    /// t_i exponents
    std::vector<double> t;
};

/// Powers of T reduced
///
/// @tparam T The basic data type
///
/// \f$ v = \displaystyle\sum_{i=0}^{n} a_i T_{r}^{t_i}\f$
template <typename TYPE>
class PowersOfTreduced {
public:
    ///
    ///
    /// @param T_critical Critical temperature [K]
    /// @param a Array of \f$a_i\f$ coefficients
    /// @param t Array of \f$t_i\f$ exponents
    PowersOfTreduced(double T_critical,
                     const std::vector<double> & a,
                     const std::vector<double> & t) :
        T_crit(T_critical),
        a(a),
        t(t)
    {
    }

    /// Evaluate the model
    ///
    /// @param T Temperature
    /// @return Computed value
    TYPE
    value(TYPE T) const
    {
        TYPE Tr = T / this->T_crit;
        TYPE sum = 0;
        for (std::size_t i = 0; i < a.size(); ++i)
            sum += this->a[i] * math::pow(Tr, this->t[i]);
        return sum;
    }

private:
    /// Critical temperature [K]
    double T_crit;
    /// a_i coefficients
    std::vector<double> a;
    /// t_i exponents
    std::vector<double> t;
};

} // namespace fprops
