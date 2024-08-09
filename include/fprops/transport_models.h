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

template <typename TYPE>
class PowersOfTr {
public:
    ///
    ///
    /// @param T_reducing Reducing temperature [K]
    /// @param a \f$a_i\f$
    /// @param t \f$t_i\f$
    PowersOfTr(double T_reducing, const std::vector<double> & a, const std::vector<double> & t) :
        T_reducing(T_reducing),
        a(a),
        t(t)
    {
    }

    TYPE
    value(TYPE T) const
    {
        TYPE Tr = T / this->T_reducing;
        TYPE summer = 0;
        for (std::size_t i = 0; i < a.size(); ++i)
            summer += this->a[i] * pow(Tr, this->t[i]);
        return summer;
    }

private:
    double T_reducing;
    std::vector<double> a;
    std::vector<double> t;
};

///
template <typename TYPE>
class FrictionTheory {
public:
    FrictionTheory(double c1,
                   double c2,
                   const std::vector<double> & Ai,
                   const std::vector<double> & Aa,
                   double Na,
                   const std::vector<double> & Ar,
                   double Nr,
                   const std::vector<double> & Aaa,
                   double Naa,
                   const std::vector<double> & Arr,
                   const std::vector<double> & Adrdr,
                   double Nrr,
                   const std::vector<double> & Aii,
                   double Nii,
                   const std::vector<double> & Arrr,
                   double Nrrr,
                   const std::vector<double> & Aaaa,
                   double Naaa) :
        c1(c1),
        c2(c2),
        Ai(Ai),
        Aa(Aa),
        Na(Na),
        Ar(Ar),
        Nr(Nr),
        Aaa(Aaa),
        Naa(Naa),
        Arr(Arr),
        Adrdr(Adrdr),
        Nrr(Nrr),
        Aii(Aii),
        Nii(Nii),
        Arrr(Arrr),
        Nrrr(Nrrr),
        Aaaa(Aaaa),
        Naaa(Naaa)
    {
    }

    TYPE
    value(double tau, double p, double pr, double pid) const
    {
        TYPE kii = 0, krrr = 0, kaaa = 0, krr, kdrdr;

        TYPE psi1 = std::exp(tau) - this->c1;
        TYPE psi2 = std::exp(math::pow(tau, 2)) - this->c2;

        TYPE ki = (this->Ai[0] + this->Ai[1] * psi1 + this->Ai[2] * psi2) * tau;

        TYPE ka =
            (this->Aa[0] + this->Aa[1] * psi1 + this->Aa[2] * psi2) * math::pow(tau, this->Na);
        TYPE kr =
            (this->Ar[0] + this->Ar[1] * psi1 + this->Ar[2] * psi2) * math::pow(tau, this->Nr);
        TYPE kaa =
            (this->Aaa[0] + this->Aaa[1] * psi1 + this->Aaa[2] * psi2) * math::pow(tau, this->Naa);
        if (this->Arr.empty()) {
            krr = 0;
            kdrdr = (this->Adrdr[0] + this->Adrdr[1] * psi1 + this->Adrdr[2] * psi2) *
                    math::pow(tau, this->Nrr);
        }
        else {
            krr = (this->Arr[0] + this->Arr[1] * psi1 + this->Arr[2] * psi2) *
                  math::pow(tau, this->Nrr);
            kdrdr = 0;
        }
        if (!this->Aii.empty()) {
            kii = (this->Aii[0] + this->Aii[1] * psi1 + this->Aii[2] * psi2) *
                  math::pow(tau, this->Nii);
        }
        if (!this->Arrr.empty() && !this->Aaaa.empty()) {
            krrr = (this->Arrr[0] + this->Arrr[1] * psi1 + this->Arrr[2] * psi2) *
                   math::pow(tau, this->Nrrr);
            kaaa = (this->Aaaa[0] + this->Aaaa[1] * psi1 + this->Aaaa[2] * psi2) *
                   math::pow(tau, this->Naaa);
        }

        TYPE pa = p - pr;
        TYPE deltapr = pr - pid;
        TYPE eta_f = ka * pa + kr * deltapr + ki * pid + kaa * pa * pa + kdrdr * deltapr * deltapr +
                     krr * pr * pr + kii * pid * pid + krrr * pr * pr * pr + kaaa * pa * pa * pa;
        return eta_f;
    }

private:
    double c1;
    double c2;

    std::vector<double> Ai;

    std::vector<double> Aa;
    double Na;

    std::vector<double> Ar;
    double Nr;

    std::vector<double> Aaa;
    double Naa;

    std::vector<double> Arr;
    std::vector<double> Adrdr;
    double Nrr;

    std::vector<double> Aii;
    double Nii;

    std::vector<double> Arrr;
    double Nrrr;

    std::vector<double> Aaaa;
    double Naaa;
};

} // namespace fprops
