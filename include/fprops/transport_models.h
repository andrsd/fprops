// SPDX-FileCopyrightText: 2022 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#pragma once

#include "fprops/numerics.h"
#include <vector>
#include <array>
#include <cassert>

namespace fprops {

/// Polynomial model
///
/// @tparam TYPE The basic data type
///
/// \f$ lambda_0 = \displaystyle\sum_{i=1}^{n} B_i \cdot \tau^{t_i} \cdot \delta^{d_i} \f$
template <typename TYPE>
class Polynomial {
public:
    /// Polynomial model
    ///
    /// @param B Coefficients B
    /// @param t Exponents t
    /// @param d Exponents d
    Polynomial(const std::vector<double> & B,
               const std::vector<double> & t,
               const std::vector<double> & d) :
        B(B),
        t(t),
        d(d)
    {
        assert(B.size() == t.size());
        assert(B.size() == d.size());
    }

    /// Evaluate the model
    ///
    /// @param delta \f$delta\f$
    /// @param tau \f$\tau\f$
    /// @return The computed value
    TYPE
    value(double tau, double delta) const
    {
        TYPE sum = 0.;
        for (std::size_t i = 0; i < this->B.size(); i++)
            sum += this->B[i] * math::pow(tau, this->t[i]) * math::pow(delta, this->d[i]);
        return sum;
    }

private:
    std::vector<double> B;
    std::vector<double> t;
    std::vector<double> d;
};

/// Polynomial-exponential model
///
/// @tparam TYPE The basic data type
///
/// \f$ lambda_0 = \displaystyle\sum_{i=1}^{n} A_i \cdot \tau^{t_i} \cdot \delta^{d_i} \cdot
/// exp(-\gamma_i*\delta^{l_i})\f$
template <typename TYPE>
class PolynomialAndExponential {
public:
    /// Polynomial model
    ///
    /// @param A Coefficients A
    /// @param t Exponents t
    /// @param d Exponents d
    /// @param gamma Coefficents gamma
    /// @param l Exponents l
    PolynomialAndExponential(const std::vector<double> & A,
                             const std::vector<double> & t,
                             const std::vector<double> & d,
                             const std::vector<double> & gamma,
                             const std::vector<double> & l) :
        A(A),
        t(t),
        d(d),
        gamma(gamma),
        l(l)
    {
        assert(A.size() == t.size());
        assert(A.size() == d.size());
        assert(A.size() == gamma.size());
        assert(A.size() == l.size());
    }

    /// Evaluate the model
    ///
    /// @param delta \f$delta\f$
    /// @param tau \f$\tau\f$
    /// @return The computed value
    TYPE
    value(double tau, double delta) const
    {
        TYPE sum = 0.;
        for (std::size_t i = 0; i < this->A.size(); ++i)
            sum += this->A[i] * math::pow(tau, this->t[i]) * math::pow(delta, this->d[i]) *
                   std::exp(-this->gamma[i] * math::pow(delta, this->l[i]));
        return sum;
    }

private:
    std::vector<double> A;
    std::vector<double> t;
    std::vector<double> d;
    std::vector<double> gamma;
    std::vector<double> l;
};

/// Model for computing viscosity (eta_0)
///
/// @tparam TYPE The basic data type
///
/// \f$ lambda_0 = A_0 * \eta_0 + \displaystyle\sum_{i=1}^{n} A_i \cdot \tau^{t_i} \f$
template <typename TYPE>
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
    TYPE
    value(double eta0, double tau) const
    {
        TYPE sum = this->A[0] * eta0 * 1e6;
        for (unsigned int i = 1; i < A.size(); i++)
            sum += this->A[i] * math::pow(tau, this->t[i]);
        return sum;
    }

private:
    std::vector<double> A;
    std::vector<double> t;
};

/// Polynomial ratio model
///
/// \f$ lambda = \frac{\displaystyle\sum_{i} A_i \cdot T_r^{n_i}}{\displaystyle\sum_{j} B_j \cdot
/// T_r^{m_j}}\f$
template <typename TYPE>
class PolynomialRatio {
public:
    PolynomialRatio(const std::vector<double> & A,
                    const std::vector<double> & n,
                    const std::vector<double> & B,
                    const std::vector<double> & m) :
        A(A),
        n(n),
        B(B),
        m(m)
    {
        assert(A.size() == n.size());
        assert(B.size() == m.size());
    }

    /// Evaluate model
    ///
    /// @param Tr Reduced temperature [K]
    /// @return The computed value
    TYPE
    value(double Tr) const
    {
        TYPE numerator = 0.;
        for (std::size_t i = 0; i < this->A.size(); i++)
            numerator += this->A[i] * math::pow(Tr, this->n[i]);
        TYPE denominator = 0.;
        for (std::size_t i = 0; i < this->B.size(); i++)
            denominator += this->B[i] * math::pow(Tr, this->m[i]);
        return numerator / denominator;
    }

private:
    std::vector<double> A;
    std::vector<double> n;
    std::vector<double> B;
    std::vector<double> m;
};

/// Lennard-Jones model for computing viscosity
///
/// @tparam TYPE The basic data type
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
template <typename TYPE>
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
    TYPE
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

/// Collision integral
template <typename TYPE>
class CollisionIntegral {
public:
    /// Collision integral
    ///
    /// @param C
    /// @param M Molar mass [kg/mol]
    /// @param epsilon_over_k [K]
    /// @param sigma_eta [m]
    CollisionIntegral(double C,
                      double M,
                      double epsilon_over_k,
                      double sigma_eta,
                      const std::vector<double> & a,
                      const std::vector<double> & t) :
        C(C),
        M(M),
        epsilon_over_k(epsilon_over_k),
        sigma_eta(sigma_eta),
        a(a),
        t(t)
    {
    }

    /// Evaluate the model
    ///
    /// @param T Temperature [K]
    /// @return Computed value
    TYPE
    value(double T) const
    {
        TYPE Tstar = T / this->epsilon_over_k;
        // 1e9 to convert from m to nm
        TYPE sigma_nm = this->sigma_eta * 1e9;
        // 1000 to convert from kg/mol to kg/kmol
        TYPE molar_mass_kgkmol = this->M * 1000;

        /// Both the collision integral \f$\mathfrak{S}^*\f$ and effective cross section
        /// \f$\Omega^{(2,2)}\f$ have the same form, in general we don't care which is used.  The
        /// are related through \f$\Omega^{(2,2)} = (5/4)\mathfrak{S}^*\f$ see Vesovic(JPCRD, 1990)
        /// for CO\f$_2\f$ for further information
        TYPE lnTstar = std::log(Tstar);
        TYPE S = 0;
        for (std::size_t i = 0; i < this->a.size(); i++)
            S += this->a[i] * math::pow(lnTstar, this->t[i]);
        S = std::exp(S);

        // The dilute gas component [Pa-s]
        return this->C * std::sqrt(molar_mass_kgkmol * T) / (math::pow<2>(sigma_nm) * S);
    }

private:
    double C;
    /// Molar mass [kg/mol]
    double M;
    double epsilon_over_k;
    double sigma_eta;
    /// Coefficients a_i
    std::vector<double> a;
    /// Exponents t_i
    std::vector<double> t;
};

/// Modified Batshinski-Hildebrand model
///
/// @tparam TYPE The basic data type
template <typename TYPE>
class ModifiedBatshinskiHildebrand {
public:
    /// Modified Batshinki-Hildebrand
    ModifiedBatshinskiHildebrand(const std::vector<double> & a,
                                 const std::vector<double> & d1,
                                 const std::vector<double> & t1,
                                 const std::vector<double> & gamma,
                                 const std::vector<double> & l,
                                 const std::vector<double> & f,
                                 const std::vector<double> & d2,
                                 const std::vector<double> & t2,
                                 const std::vector<double> & g,
                                 const std::vector<double> & h,
                                 const std::vector<double> & p,
                                 const std::vector<double> & q) :
        a(a),
        d1(d1),
        t1(t1),
        gamma(gamma),
        l(l),
        f(f),
        d2(d2),
        t2(t2),
        g(g),
        h(h),
        p(p),
        q(q)
    {
        assert(this->a.size() == this->d1.size());
        assert(this->a.size() == this->t1.size());
        assert(this->a.size() == this->gamma.size());
        assert(this->a.size() == this->l.size());
        assert(this->f.size() == this->d2.size());
        assert(this->f.size() == this->t2.size());
        assert(this->g.size() == this->h.size());
        assert(this->p.size() == this->q.size());
    }

    /// Evaluate the model
    ///
    /// @param delta \f$\delta\f$
    /// @param tau \f$\tau\f$
    /// @return Computed value
    TYPE
    value(double delta, double tau) const
    {
        // The first term that is formed of powers of tau (Tc/T) and delta (rho/rhoc)
        TYPE S = 0;
        for (std::size_t i = 0; i < this->a.size(); ++i)
            S += this->a[i] * math::pow(delta, this->d1[i]) * math::pow(tau, this->t1[i]) *
                 std::exp(this->gamma[i] * math::pow(delta, this->l[i]));

        // For the terms that multiplies the bracketed term with delta and delta0
        TYPE F = 0;
        for (std::size_t i = 0; i < this->f.size(); ++i)
            F += this->f[i] * math::pow(delta, this->d2[i]) * math::pow(tau, this->t2[i]);

        // for delta_0
        TYPE numer = 0;
        for (std::size_t i = 0; i < this->g.size(); i++)
            numer += this->g[i] * pow(tau, this->h[i]);
        TYPE denom = 0;
        for (std::size_t i = 0; i < this->p.size(); i++)
            denom += this->p[i] * pow(tau, this->q[i]);
        TYPE delta0 = numer / denom;

        // The higher-order-term component [Pa-s]
        return S + F * (1. / (delta0 - delta) - 1 / delta0);
    }

private:
    std::vector<double> a;
    std::vector<double> d1;
    std::vector<double> t1;
    std::vector<double> gamma;
    std::vector<double> l;
    std::vector<double> f;
    std::vector<double> d2;
    std::vector<double> t2;
    std::vector<double> g;
    std::vector<double> h;
    std::vector<double> p;
    std::vector<double> q;
};

/// Powers of temperature
///
/// @tparam TYPE The basic data type
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

/// Rainwater-Friend transport model
///
template <typename TYPE>
class RainwaterFriend {
public:
    /// Rainwater-Friend transport model
    ///
    /// @param epsilon_over_k \f$\epsilon\over k\f$ \f$[K]\f$
    /// @param sigma_eta \f$\sigma \eta\f$ \f$[m]\f$
    /// @param b Coefficients \f$b_i\f$
    /// @param t Exponents \f$t_i\f$
    RainwaterFriend(double epsilon_over_k,
                    double sigma_eta,
                    const std::vector<double> & b,
                    const std::vector<double> & t) :
        epsilon_over_k(epsilon_over_k),
        sigma_eta(sigma_eta),
        b(b),
        t(t)
    {
    }

    /// Evaluate the model
    ///
    /// @param T Temperature [K]
    /// @return Computed value [m^3/mol]
    TYPE
    value(TYPE T) const
    {
        TYPE Tstar = T / this->epsilon_over_k;
        TYPE sigma = this->sigma_eta;
        // [no units]
        TYPE B_eta_star = 0;
        for (std::size_t i = 0; i < this->b.size(); ++i)
            B_eta_star += b[i] * math::pow(Tstar, this->t[i]);
        auto B_eta = 6.02214129e23 * math::pow<3>(sigma) * B_eta_star;
        return B_eta;
    }

private:
    /// Units are [K]
    double epsilon_over_k;
    /// Units are [m]
    double sigma_eta;
    /// Coefficients b_i
    std::vector<double> b;
    /// Exponents t_i
    std::vector<double> t;
};

/// Quadratic general friction theory model
template <typename TYPE>
class QuadraticGeneralFrictionTheory {
public:
    /// Quadratic general friction theory (FT) model
    ///
    /// @param a Coefficients a [(Pa*s)/bar]
    /// @param b Coefficients b [(Pa*s)/bar]
    /// @param c Coefficients c [(Pa*s)/bar]
    /// @param A Coefficients A [(Pa*s)/bar^2]
    /// @param B Coefficients B [(Pa*s)/bar^2]
    /// @param C Coefficients C [(Pa*s)/bar^2]
    QuadraticGeneralFrictionTheory(const std::array<double, 3> & a,
                                   const std::array<double, 3> & b,
                                   const std::array<double, 3> & c,
                                   const std::array<double, 3> & A,
                                   const std::array<double, 3> & B,
                                   const std::array<double, 3> & C) :
        a(a),
        b(b),
        c(c),
        A(A),
        B(B),
        C(C)
    {
    }

    /// Evaluate the model
    ///
    /// @param tau
    /// @param p_a Pressure [bar]
    /// @param p_r Pressure [bar]
    /// @param p_id Pressure [bar]
    TYPE
    value(double tau, double p_a, double p_r, double p_id) const
    {
        // to match the convention in the paper
        TYPE Gamma = tau;
        TYPE psi1 = std::exp(tau) - 1;
        TYPE psi2 = std::exp(math::pow(tau, 2)) - 1;
        std::array<TYPE, 3> psi = { 1., psi1, psi2 };
        TYPE kappa_a = dot(psi, this->a) * Gamma;
        TYPE kappa_aa = dot(psi, this->A) * math::pow<3>(Gamma);
        TYPE kappa_r = dot(psi, this->b) * Gamma;
        TYPE kappa_rr = dot(psi, this->B) * math::pow<3>(Gamma);
        TYPE kappa_i = dot(psi, this->c) * Gamma;
        TYPE kappa_ii = dot(psi, this->C) * math::pow<3>(Gamma);
        TYPE Delta_p_r = p_r - p_id;
        return kappa_i * p_id + kappa_r * Delta_p_r + kappa_a * p_a +
               kappa_ii * math::pow<2>(p_id) + kappa_rr * math::pow<2>(Delta_p_r) +
               kappa_aa * math::pow<2>(p_a);
    }

private:
    template <auto N>
    inline TYPE
    dot(const std::array<TYPE, N> & a, const std::array<TYPE, N> & b) const
    {
        TYPE prod = 0.;
        for (auto i = 0; i < N; i++)
            prod += a[i] * b[i];
        return prod;
    }

    std::array<double, 3> a;
    std::array<double, 3> b;
    std::array<double, 3> c;
    std::array<double, 3> A;
    std::array<double, 3> B;
    std::array<double, 3> C;
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
