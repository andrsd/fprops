// SPDX-FileCopyrightText: 2022 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#pragma once

#include "fprops/single_phase_fluid_properties.h"
#include "fprops/numerics.h"
#include <cmath>
#include <vector>
#include <cassert>

namespace fprops {

/// Base class for fluid properties based on Helmholtz equation of state
///
/// This class is based on `HelmholtzFluidProperties.h` from `idaholab/moose/fluid_properties`
/// module
class Helmholtz : public SinglePhaseFluidProperties {
public:
    /// @param R Universal gas constant \f$[J/(mol-K)]\f$
    /// @param M Molar mass \f$[kg/mol]\f$
    /// @param rho_c Critical density \f$[kg/m^3]\f$
    /// @param T_c Critical temperature \f$[K]\f$
    Helmholtz(double R, double M, double rho_c, double T_c);

    [[nodiscard]] State rho_T(double rho, double T) const override;
    [[nodiscard]] State rho_p(double rho, double p) const override;
    [[nodiscard]] State p_T(double p, double T) const override;
    [[nodiscard]] State v_u(double v, double u) const override;
    [[nodiscard]] State h_s(double h, double s) const override;

private:
    /// Helmholtz free energy
    ///
    /// @param delta Scaled density \f$[-]\f$
    /// @param tau Scaled temperature \f$[-]\f$
    /// @return Helmholtz free energy (\f$\alpha\f$)
    [[nodiscard]] virtual double alpha(double delta, double tau) const = 0;

    /// Derivative of Helmholtz free energy wrt delta
    ///
    /// @param delta Scaled density \f$[-]\f$
    /// @param tau Scaled temperature \f$[-]\f$
    /// @return Derivative of Helmholtz free energy wrt delta
    [[nodiscard]] virtual double dalpha_ddelta(double delta, double tau) const = 0;

    /// Derivative of Helmholtz free energy wrt tau
    ///
    /// @param delta Scaled density \f$[-]\f$
    /// @param tau Scaled temperature \f$[-]\f$
    /// @return Derivative of Helmholtz free energy wrt tau
    [[nodiscard]] virtual double dalpha_dtau(double delta, double tau) const = 0;

    /// Second derivative of Helmholtz free energy wrt delta
    ///
    /// @param delta Scaled density \f$[-]\f$
    /// @param tau Scaled temperature \f$[-]\f$
    /// @return Second derivative of Helmholtz free energy wrt delta
    [[nodiscard]] virtual double d2alpha_ddelta2(double delta, double tau) const = 0;

    /// Second derivative of Helmholtz free energy wrt tau
    ///
    /// @param delta Scaled density \f$[-]\f$
    /// @param tau Scaled temperature \f$[-]\f$
    /// @return Second derivative of Helmholtz free energy wrt tau
    [[nodiscard]] virtual double d2alpha_dtau2(double delta, double tau) const = 0;

    /// Second derivative of Helmholtz free energy wrt delta and tau
    ///
    /// @param delta Scaled density \f$[-]\f$
    /// @param tau Scaled temperature \f$[-]\f$
    /// @return Second derivative of Helmholtz free energy wrt delta and tau
    [[nodiscard]] virtual double d2alpha_ddeltatau(double delta, double tau) const = 0;

    /// Density given pressure and temperature
    ///
    /// @param p Pressure \f$[Pa]\f$
    /// @param T Temperature \f$[K]\f$
    /// @return Density \f$[kg/m^3]\f$
    [[nodiscard]] double rho_from_p_T(double p, double T) const;

    /// Temperature given density and pressure
    ///
    /// @param rho Density \f$[kg/m^3]\f$
    /// @param p Pressure \f$[Pa]\f$
    /// @return Temperature \f$[K]\f$
    [[nodiscard]] double T_from_rho_p(double rho, double p) const;

    /// Tau given specific volume and internal energy
    ///
    /// @param v Specific volume \f$[m^3/kg]\f$
    /// @param u Specific internal energy \f$[J/kg]\f$
    /// @return Tau (scaled temperature)
    [[nodiscard]] double tau_from_v_u(double v, double u) const;

    /// Dynamic viscosity
    ///
    /// @param rho Density \f$[kg/m^3]\f$
    /// @param T Temperature \f$[K]\f$
    /// @return Dynamic viscosity \f$[Pa-s]\f$
    [[nodiscard]] virtual double mu_from_rho_T(double rho, double T) const = 0;

    /// Thermal conductivity
    ///
    /// @param rho Density \f$[kg/m^3]\f$
    /// @param T Temperature \f$[K]\f$
    /// @return Thermal conductivity \f$[W/(m-K)]\f$
    [[nodiscard]] virtual double k_from_rho_T(double rho, double T) const = 0;

    /// Universal gas constant \f$[J/(mol-K)]\f$
    const double R;
    /// Molar mass \f$[kg/mol]\f$
    const double M;
    /// Critical density \f$[kg/m^3]\f$
    const double rho_c;
    /// Critical temperature \f$[K]\f$
    const double T_c;

protected:
    /// The leading term in the EOS used to set the desired reference state
    ///
    /// @tparam T The basic data type
    ///
    /// \f$\alpha = \ln(\delta) + a_1 + a_2 \tau\f$
    template <typename T>
    class IdealGasLead {
    public:
        IdealGasLead(double a1, double a2) : a1(a1), a2(a2) {}

        T
        alpha(T delta, T tau) const
        {
            return std::log(delta) + this->a1 + this->a2 * tau;
        }

        T
        ddelta(T delta, T tau) const
        {
            return 1.0 / delta;
        }

        T
        dtau(T delta, T tau) const
        {
            return this->a2;
        }

        T
        d2delta(T delta, T tau) const
        {
            return -1.0 / delta / delta;
        }

    private:
        double a1, a2;
    };

    /// Log(tau) term
    ///
    /// @tparam T The basic data type
    ///
    /// \f$\alpha = a_1\ln(\tau)\f$
    template <typename T>
    class IdealGasLogTau {
    public:
        explicit IdealGasLogTau(double a1) : a1(a1) {}

        T
        alpha(T delta, T tau) const
        {
            return this->a1 * std::log(tau);
        }

        T
        dtau(T delta, T tau) const
        {
            return this->a1 / tau;
        }

        T
        d2tau(T delta, T tau) const
        {
            return -this->a1 / tau / tau;
        }

    private:
        double a1;
    };

    /// Sum of powers used for ideal part of \f$\alpha\f$
    ///
    /// @tparam T The basic data type
    ///
    /// \f$\alpha = \displaystyle\sum_{i=0}^n N_i \tau^{t_i}\f$
    template <typename T>
    class IdealGasPower {
    public:
        IdealGasPower(const std::vector<T> & n, const std::vector<T> & t) : n(n), t(t)
        {
            assert(n.size() == t.size());
        }

        T
        alpha(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); ++i)
                sum += this->n[i] * math::pow(tau, this->t[i]);
            return sum;
        }

        T
        dtau(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < n.size(); ++i)
                sum += this->n[i] * this->t[i] * math::pow(tau, this->t[i] - 1);
            return sum;
        }

        T
        d2tau(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); ++i)
                sum += this->n[i] * this->t[i] * (this->t[i] - 1) * math::pow(tau, this->t[i] - 2);
            return sum;
        }

    private:
        std::vector<T> n;
        std::vector<T> t;
    };

    /// Offset for enthalpy and entropy
    ///
    /// @tparam T The basic data type
    ///
    /// \f$\alpha = a_1 + a_2 \tau\f$
    template <typename T>
    class IdealEnthalpyEntropyOffset {
    public:
        IdealEnthalpyEntropyOffset(double a1, double a2) : a1(a1), a2(a2) {}

        T
        alpha(T delta, T tau) const
        {
            return this->a1 + this->a2 * tau;
        }

        T
        dtau(T delta, T tau) const
        {
            return this->a2;
        }

    private:
        double a1;
        double a2;
    };

    /// Generalized Planck-Einstein model
    ///
    /// @tparam T The basic data type
    ///
    /// \f$ \alpha = \displaystyle\sum_{i=0}^{n} N_i \log(c_i + d_i \exp(\theta \tau)) \f$
    template <typename T>
    class IdealPlanckEinsteinGeneralized {
    public:
        IdealPlanckEinsteinGeneralized(const std::vector<T> & n,
                                       const std::vector<T> & theta,
                                       const std::vector<T> & c,
                                       const std::vector<T> & d) :
            n(n),
            theta(theta),
            c(c),
            d(d)
        {
            assert(n.size() == theta.size());
            assert(n.size() == c.size());
            assert(n.size() == d.size());
        }

        T
        alpha(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum +=
                    this->n[i] * std::log(this->c[i] + this->d[i] * std::exp(this->theta[i] * tau));
            return sum;
        }

        T
        dtau(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); ++i) {
                double exp_theta_tau = std::exp(this->theta[i] * tau);
                sum += this->n[i] * this->theta[i] * this->d[i] * exp_theta_tau /
                       (this->c[i] + this->d[i] * exp_theta_tau);
            }
            return sum;
        }

        T
        d2tau(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); ++i) {
                double exp_theta_tau = std::exp(this->theta[i] * tau);
                sum += this->n[i] * math::pow<2>(this->theta[i]) * this->c[i] * this->d[i] *
                       exp_theta_tau / math::pow<2>(this->c[i] + this->d[i] * exp_theta_tau);
            }
            return sum;
        }

    private:
        std::vector<T> n;
        std::vector<T> theta;
        std::vector<T> c;
        std::vector<T> d;
    };

    /// Planck-Einstein model as a function of temperature
    ///
    /// @tparam T The basic data type
    ///
    /// \f$ \alpha = \displaystyle\sum_{i=0}^{n} N_i \log(1 - \exp({v_i \over T_{crit}} \tau)) \f$
    template <typename T>
    class IdealGasPlanckEinsteinFunctionT {
    public:
        IdealGasPlanckEinsteinFunctionT(double T_crit,
                                        const std::vector<T> & n,
                                        const std::vector<T> & v)
        {
            std::vector<T> theta(n.size(), 0.);
            for (std::size_t i = 0; i < v.size(); ++i)
                theta[i] = -v[i] / T_crit;
            std::vector<T> c(n.size(), 1);
            std::vector<T> d(c.size(), -1);
            this->generalized = new IdealPlanckEinsteinGeneralized<T>(n, theta, c, d);
        }

        ~IdealGasPlanckEinsteinFunctionT() { delete this->generalized; }

        T
        alpha(T delta, T tau) const
        {
            return this->generalized->alpha(delta, tau);
        }

        T
        dtau(T delta, T tau) const
        {
            return this->generalized->dtau(delta, tau);
        }

        T
        d2tau(T delta, T tau) const
        {
            return this->generalized->d2tau(delta, tau);
        }

    private:
        IdealPlanckEinsteinGeneralized<T> * generalized;
    };

    /// Planck-Einstein
    ///
    /// @tparam T basic data type
    ///
    /// \f$ \alpha = \displaystyle\sum_{i=0}^{n} N_i \log(1 - \exp(-t_i \tau)) \f$
    template <typename T>
    class IdealGasPlanckEinstein {
    public:
        IdealGasPlanckEinstein(const std::vector<T> & n, const std::vector<T> & t)
        {
            assert(n.size() == t.size());
            std::vector<T> theta(n.size(), 0.);
            for (std::size_t i = 0; i < t.size(); i++)
                theta[i] = -t[i];
            std::vector<T> c(n.size(), 1);
            std::vector<T> d(n.size(), -1);
            this->generalized = new IdealPlanckEinsteinGeneralized<T>(n, theta, c, d);
        }

        ~IdealGasPlanckEinstein() { delete this->generalized; }

        T
        alpha(T delta, T tau) const
        {
            return this->generalized->alpha(delta, tau);
        }

        T
        dtau(T delta, T tau) const
        {
            return this->generalized->dtau(delta, tau);
        }

        T
        d2tau(T delta, T tau) const
        {
            return this->generalized->d2tau(delta, tau);
        }

    private:
        IdealPlanckEinsteinGeneralized<T> * generalized;
    };

    /// Power model
    ///
    /// @tparam T basic data type
    ///
    /// \f$ \alpha = \displaystyle\sum_{i=0}^{n} N_i \delta^{d_i} \tau^{t_i}  \f$
    template <typename T>
    class ResidualPower {
    public:
        ResidualPower(const std::vector<T> & n,
                      const std::vector<T> & d,
                      const std::vector<T> & t) :
            n(n),
            d(d),
            t(t)
        {
            assert(n.size() == d.size());
            assert(n.size() == t.size());
        }

        T
        alpha(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum += this->n[i] * math::pow(delta, this->d[i]) * math::pow(tau, this->t[i]);
            return sum;
        }

        T
        ddelta(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum += this->n[i] * this->d[i] * math::pow(delta, this->d[i] - 1) *
                       math::pow(tau, this->t[i]);
            return sum;
        }

        T
        dtau(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum += this->n[i] * math::pow(delta, this->d[i]) * this->t[i] *
                       math::pow(tau, this->t[i] - 1);
            return sum;
        }

        T
        d2delta(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum += this->n[i] * this->d[i] * (this->d[i] - 1) *
                       math::pow(delta, this->d[i] - 2) * math::pow(tau, this->t[i]);
            return sum;
        }

        T
        d2tau(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum += this->n[i] * math::pow(delta, this->d[i]) * this->t[i] * (this->t[i] - 1) *
                       math::pow(tau, this->t[i] - 2);
            return sum;
        }

        T
        d2deltatau(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum += this->n[i] * this->d[i] * math::pow(delta, this->d[i] - 1) * this->t[i] *
                       math::pow(tau, this->t[i] - 1);
            return sum;
        }

    private:
        std::vector<T> n;
        std::vector<T> d;
        std::vector<T> t;
    };

    /// Power model with exponential term
    ///
    /// @tparam T basic data type
    /// @tparam L type of the `l`-coefficient
    ///
    /// \f$ \alpha = \displaystyle\sum_{i=0}^{n} N_i \delta^{d_i} \tau^{t_i} \exp(-\delta^{l_i}) \f$
    template <typename T, typename L>
    class ResidualPowerExp {
    public:
        ResidualPowerExp(const std::vector<T> & n,
                         const std::vector<T> & d,
                         const std::vector<T> & t,
                         const std::vector<L> & l) :
            n(n),
            d(d),
            t(t),
            l(l)
        {
            assert(n.size() == d.size());
            assert(n.size() == t.size());
            assert(n.size() == l.size());
        }

        T
        alpha(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum += this->n[i] * std::pow(delta, this->d[i]) * std::pow(tau, this->t[i]) *
                       std::exp(-std::pow(delta, this->l[i]));
            return sum;
        }

        T
        ddelta(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum += this->n[i] * std::pow(delta, this->d[i] - 1) * std::pow(tau, this->t[i]) *
                       std::exp(-std::pow(delta, this->l[i])) *
                       (this->d[i] - this->l[i] * std::pow(delta, this->l[i]));
            return sum;
        }

        T
        dtau(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum += this->n[i] * this->t[i] * std::pow(delta, this->d[i]) *
                       std::pow(tau, this->t[i] - 1) * std::exp(-std::pow(delta, this->l[i]));
            return sum;
        }

        T
        d2delta(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum +=
                    this->n[i] * std::pow(delta, this->d[i] - 2) * std::pow(tau, this->t[i]) *
                    std::exp(-std::pow(delta, this->l[i])) *
                    (math::pow<2>(this->l[i]) * std::pow(delta, 2 * this->l[i]) +
                     (this->d[i] - 1) * this->d[i] -
                     this->l[i] * (2 * this->d[i] + this->l[i] - 1) * std::pow(delta, this->l[i]));
            return sum;
        }

        T
        d2tau(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum += this->n[i] * (this->t[i] - 1) * this->t[i] * std::pow(delta, this->d[i]) *
                       std::exp(-std::pow(delta, this->l[i])) * std::pow(tau, this->t[i] - 2);
            return sum;
        }

        T
        d2deltatau(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum += this->n[i] * this->t[i] * std::pow(delta, this->d[i] - 1) *
                       std::pow(tau, this->t[i] - 1) * std::exp(-std::pow(delta, this->l[i])) *
                       (this->d[i] - this->l[i] * std::pow(delta, this->l[i]));
            return sum;
        }

    private:
        std::vector<T> n;
        std::vector<T> d;
        std::vector<T> t;
        std::vector<L> l;
    };

    /// Gaussian model for residual part
    ///
    /// @tparam T The basic data type
    ///
    /// \f$ \alpha = \displaystyle\sum_{i=0}^n N_i \delta^{d_i} \tau^{t_i}
    /// \exp(-\eta_i (\delta - \epsilon_i)^2 - \beta_i (\tau - \gamma_i)^2) \f$
    template <typename T>
    class ResidualGaussian {
    public:
        ResidualGaussian(const std::vector<T> & n,
                         const std::vector<T> & d,
                         const std::vector<T> & t,
                         const std::vector<T> & eta,
                         const std::vector<T> & epsilon,
                         const std::vector<T> & beta,
                         const std::vector<T> & gamma) :
            n(n),
            d(d),
            t(t),
            eta(eta),
            epsilon(epsilon),
            beta(beta),
            gamma(gamma)
        {
            assert(n.size() == d.size());
            assert(n.size() == t.size());
            assert(n.size() == eta.size());
            assert(n.size() == epsilon.size());
            assert(n.size() == beta.size());
            assert(n.size() == gamma.size());
        }

        T
        alpha(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum += this->n[i] * std::pow(delta, this->d[i]) * std::pow(tau, this->t[i]) *
                       std::exp(-this->eta[i] * math::pow<2>(delta - this->epsilon[i]) -
                                this->beta[i] * math::pow<2>(tau - this->gamma[i]));
            return sum;
        }

        T
        ddelta(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum += this->n[i] * std::pow(delta, this->d[i] - 1) * std::pow(tau, this->t[i]) *
                       (this->d[i] - 2.0 * delta * this->eta[i] * (delta - this->epsilon[i])) *
                       std::exp(-this->eta[i] * math::pow<2>(delta - this->epsilon[i]) -
                                this->beta[i] * math::pow<2>(tau - this->gamma[i]));
            return sum;
        }

        T
        dtau(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum += this->n[i] * std::pow(delta, this->d[i]) * std::pow(tau, this->t[i] - 1) *
                       (2.0 * this->beta[i] * (this->gamma[i] - tau) * tau + this->t[i]) *
                       std::exp(-this->eta[i] * math::pow<2>(delta - this->epsilon[i]) -
                                this->beta[i] * math::pow<2>(tau - this->gamma[i]));
            return sum;
        }

        T
        d2delta(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum += this->n[i] * std::pow(delta, this->d[i] - 2) * std::pow(tau, this->t[i]) *
                       (2 * math::pow<2>(delta) * this->eta[i] *
                            (2 * this->eta[i] * math::pow<2>(delta - this->epsilon[i]) - 1) +
                        math::pow<2>(this->d[i]) +
                        this->d[i] * (4 * delta * this->eta[i] * (this->epsilon[i] - delta) - 1)) *
                       std::exp(-this->eta[i] * math::pow<2>(delta - this->epsilon[i]) -
                                this->beta[i] * math::pow<2>(tau - this->gamma[i]));
            return sum;
        }

        T
        d2tau(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum += this->n[i] * std::pow(delta, this->d[i]) * std::pow(tau, this->t[i] - 2) *
                       (2 * this->beta[i] * math::pow<2>(tau) *
                            (2 * this->beta[i] * math::pow<2>(this->gamma[i] - tau) - 1) +
                        math::pow<2>(this->t[i]) +
                        this->t[i] * (4 * this->beta[i] * (this->gamma[i] - tau) * tau - 1)) *
                       std::exp(-this->eta[i] * math::pow<2>(delta - this->epsilon[i]) -
                                this->beta[i] * math::pow<2>(tau - this->gamma[i]));
            return sum;
        }

        T
        d2deltatau(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum += this->n[i] * std::pow(delta, this->d[i] - 1) *
                       std::pow(tau, this->t[i] - 1) *
                       (2. * this->beta[i] * (this->gamma[i] - tau) * tau + this->t[i]) *
                       (2. * delta * this->eta[i] * (this->epsilon[i] - delta) + this->d[i]) *
                       std::exp(-this->beta[i] * math::pow<2>(this->gamma[i] - tau) -
                                this->eta[i] * math::pow<2>(delta - this->epsilon[i]));
            return sum;
        }

    private:
        std::vector<T> n;
        std::vector<T> d;
        std::vector<T> t;
        std::vector<T> eta;
        std::vector<T> epsilon;
        std::vector<T> beta;
        std::vector<T> gamma;
    };

    template <typename T>
    class ResidualNonAnalytic {
    public:
        ResidualNonAnalytic(const std::vector<T> & n,
                            const std::vector<T> & a,
                            const std::vector<T> & b,
                            const std::vector<T> & beta,
                            const std::vector<T> & A,
                            const std::vector<T> & B,
                            const std::vector<T> & C,
                            const std::vector<T> & D) :
            n(n),
            a(a),
            b(b),
            beta(beta),
            A(A),
            B(B),
            C(C),
            D(D)
        {
            assert(n.size() == a.size());
            assert(n.size() == b.size());
            assert(n.size() == beta.size());
            assert(n.size() == A.size());
            assert(n.size() == B.size());
            assert(n.size() == C.size());
            assert(n.size() == D.size());
        }

        T
        alpha(T delta, T tau) const
        {
            fiddle(delta, tau);
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum += delta * this->n[i] * DELTA_bi(i, delta, tau) * PSI(i, delta, tau);
            return sum;
        }

        T
        ddelta(T delta, T tau) const
        {
            fiddle(delta, tau);
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++) {
                sum += this->n[i] * (DELTA_bi(i, delta, tau) * PSI(i, delta, tau) +
                                     delta * dDELTA_bi_ddelta(i, delta, tau) * PSI(i, delta, tau) +
                                     delta * DELTA_bi(i, delta, tau) * dPSI_ddelta(i, delta, tau));
            }
            return sum;
        }

        T
        dtau(T delta, T tau) const
        {
            fiddle(delta, tau);
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++) {
                sum += this->n[i] * delta *
                       (dDELTA_bi_dtau(i, delta, tau) * PSI(i, delta, tau) +
                        DELTA_bi(i, delta, tau) * dPSI_dtau(i, delta, tau));
            }
            return sum;
        }

        T
        d2delta(T delta, T tau) const
        {
            fiddle(delta, tau);
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++) {
                sum += this->n[i] *
                       (dDELTA_bi_ddelta(i, delta, tau) * PSI(i, delta, tau) +
                        DELTA_bi(i, delta, tau) * dPSI_ddelta(i, delta, tau) +
                        dDELTA_bi_ddelta(i, delta, tau) * PSI(i, delta, tau) +
                        delta * d2DELTA_bi_ddelta2(i, delta, tau) * PSI(i, delta, tau) +
                        delta * dDELTA_bi_ddelta(i, delta, tau) * dPSI_ddelta(i, delta, tau) +
                        DELTA_bi(i, delta, tau) * dPSI_ddelta(i, delta, tau) +
                        delta * dDELTA_bi_ddelta(i, delta, tau) * dPSI_ddelta(i, delta, tau) +
                        delta * DELTA_bi(i, delta, tau) * d2PSI_ddelta2(i, delta, tau));
            }
            return sum;
        }

        T
        d2tau(T delta, T tau) const
        {
            fiddle(delta, tau);
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++) {
                sum += this->n[i] * delta *
                       (d2DELTA_bi_dtau2(i, delta, tau) * PSI(i, delta, tau) +
                        2 * dDELTA_bi_dtau(i, delta, tau) * dPSI_dtau(i, delta, tau) +
                        DELTA_bi(i, delta, tau) * d2PSI_dtau2(i, delta, tau));
            }
            return sum;
        }

        T
        d2deltatau(T delta, T tau) const
        {
            fiddle(delta, tau);
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++) {
                sum += this->n[i] *
                       (dDELTA_bi_dtau(i, delta, tau) * PSI(i, delta, tau) +
                        delta * (d2DELTA_bi_ddeltatau(i, delta, tau) * PSI(i, delta, tau) +
                                 dDELTA_bi_dtau(i, delta, tau) * dPSI_ddelta(i, delta, tau) +
                                 dDELTA_bi_ddelta(i, delta, tau) * dPSI_dtau(i, delta, tau) +
                                 DELTA_bi(i, delta, tau) * d2PSI_ddeltatau(i, delta, tau)));
            }
            return sum;
        }

    private:
        // Manipulate `delta` and `tau` such that we are not evaluating at the critical point but
        // very close to it
        void
        fiddle(T & delta, T & tau) const
        {
            T EPSILON = std::numeric_limits<T>::epsilon();
            if (std::abs(tau - 1) < 10 * EPSILON)
                tau = 1.0 + 10 * EPSILON;
            if (std::abs(delta - 1) < 10 * EPSILON)
                delta = 1.0 + 10 * EPSILON;
        }

        T
        theta(std::size_t i, T delta, T tau) const
        {
            return (1.0 - tau) +
                   this->A[i] * std::pow(math::pow<2>(delta - 1.0), 1.0 / (2.0 * this->beta[i]));
        }

        T
        dtheta_ddelta(std::size_t i, T delta, T tau) const
        {
            return this->A[i] * std::pow(math::pow<2>(delta - 1), 1. / (2. * this->beta[i])) /
                   (this->beta[i] * (delta - 1));
        }

        T
        dtheta_dtau(std::size_t i, T delta, T tau) const
        {
            return -1;
        }

        T
        d2theta_ddelta2(std::size_t i, T delta, T tau) const
        {
            return -1 *
                   (this->A[i] * (this->beta[i] - 1) *
                    std::pow(math::pow<2>(delta - 1), 1. / (2 * this->beta[i]) - 1.)) /
                   (math::pow<2>(this->beta[i]));
        }

        T
        DELTA(std::size_t i, T delta, T tau) const
        {
            return math::pow<2>(theta(i, delta, tau)) +
                   this->B[i] * std::pow(math::pow<2>(delta - 1.0), this->a[i]);
        }

        T
        dDELTA_ddelta(std::size_t i, T delta, T tau) const
        {
            return 2. * theta(i, delta, tau) * dtheta_ddelta(i, delta, tau) +
                   (2 * this->a[i] * this->B[i] * std::pow(delta - 1, 2 * this->a[i] - 1));
        }

        T
        dDELTA_dtau(std::size_t i, T delta, T tau) const
        {
            return 2. * theta(i, delta, tau) * dtheta_dtau(i, delta, tau);
        }

        T
        d2DELTA_ddelta2(std::size_t i, T delta, T tau) const
        {
            return 2. * (math::pow<2>(dtheta_ddelta(i, delta, tau) +
                                      theta(i, delta, tau) * d2theta_ddelta2(i, delta, tau))) +
                   2 * this->a[i] * (2 * this->a[i] - 1) * this->B[i] *
                       std::pow(delta - 1, 2 * this->a[i] - 2);
        }

        T
        d2DELTA_dtau2(std::size_t i, T delta, T tau) const
        {
            return 2. * math::pow<2>(dtheta_dtau(i, delta, tau));
        }

        T
        d2DELTA_ddeltatau(std::size_t i, T delta, T tau) const
        {
            return 2 * dtheta_ddelta(i, delta, tau) * dtheta_dtau(i, delta, tau);
        }

        T
        DELTA_bi(std::size_t i, T delta, T tau) const
        {
            return std::pow(DELTA(i, delta, tau), this->b[i]);
        }

        T
        dDELTA_bi_ddelta(std::size_t i, T delta, T tau) const
        {
            return this->b[i] * std::pow(DELTA(i, delta, tau), this->b[i] - 1) *
                   dDELTA_ddelta(i, delta, tau);
        }

        T
        dDELTA_bi_dtau(std::size_t i, T delta, T tau) const
        {
            return this->b[i] * std::pow(DELTA(i, delta, tau), this->b[i] - 1) *
                   dDELTA_dtau(i, delta, tau);
        }

        T
        d2DELTA_bi_ddelta2(std::size_t i, T delta, T tau) const
        {
            return this->b[i] * std::pow(DELTA(i, delta, tau), this->b[i] - 1) *
                   d2DELTA_ddelta2(i, delta, tau);
        }

        T
        d2DELTA_bi_dtau2(std::size_t i, T delta, T tau) const
        {
            return this->b[i] * std::pow(DELTA(i, delta, tau), this->b[i] - 1) *
                   d2DELTA_dtau2(i, delta, tau);
        }

        T
        d2DELTA_bi_ddeltatau(std::size_t i, T delta, T tau) const
        {
            return this->b[i] * std::pow(DELTA(i, delta, tau), this->b[i] - 1) *
                   d2DELTA_ddeltatau(i, delta, tau);
        }

        T
        PSI(std::size_t i, T delta, T tau) const
        {
            return std::exp(-this->C[i] * math::pow<2>(delta - 1.0) -
                            this->D[i] * math::pow<2>(tau - 1.0));
        }

        T
        dPSI_ddelta(std::size_t i, T delta, T tau) const
        {
            return -2 * this->C[i] * (delta - 1) *
                   std::exp(-this->C[i] * math::pow<2>(delta - 1) -
                            this->D[i] * math::pow<2>(tau - 1));
        }

        T
        dPSI_dtau(std::size_t i, T delta, T tau) const
        {
            return -2 * this->D[i] * (tau - 1) *
                   std::exp(-this->C[i] * math::pow<2>(delta - 1) -
                            this->D[i] * math::pow<2>(tau - 1));
        }

        T
        d2PSI_ddelta2(std::size_t i, T delta, T tau) const
        {
            return 2 * this->C[i] * (2 * this->C[i] * math::pow<2>(delta - 1) - 1) *
                   std::exp(-this->C[i] * math::pow<2>(delta - 1) -
                            this->D[i] * math::pow<2>(tau - 1));
        }

        T
        d2PSI_dtau2(std::size_t i, T delta, T tau) const
        {
            return 2 * this->D[i] * (2 * this->D[i] * math::pow<2>(tau - 1) - 1) *
                   std::exp(-this->C[i] * math::pow<2>(delta - 1) -
                            this->D[i] * math::pow<2>(tau - 1));
        }

        T
        d2PSI_ddeltatau(std::size_t i, T delta, T tau) const
        {
            return 4. * this->C[i] * (delta - 1) * this->D[i] * (tau - 1) *
                   std::exp(-this->C[i] * math::pow<2>(delta - 1) -
                            this->D[i] * math::pow<2>(tau - 1));
        }

        std::vector<T> n;
        std::vector<T> a;
        std::vector<T> b;
        std::vector<T> beta;
        std::vector<T> A;
        std::vector<T> B;
        std::vector<T> C;
        std::vector<T> D;
    };

    template <typename TYPE>
    class ResidualGaoB {
    public:
        ResidualGaoB(const std::vector<double> & n,
                     const std::vector<double> & d,
                     const std::vector<double> & t,
                     const std::vector<double> & eta,
                     const std::vector<double> & beta,
                     const std::vector<double> & gamma,
                     const std::vector<double> & epsilon,
                     const std::vector<double> & b) :
            n(n),
            d(d),
            t(t),
            eta(eta),
            beta(beta),
            gamma(gamma),
            epsilon(epsilon),
            b(b)
        {
        }

        TYPE
        alpha(TYPE delta, TYPE tau) const
        {
            TYPE a = 0.;
            for (std::size_t i = 0; i < this->n.size(); i++)
                a += this->n[i] * f_tau(i, tau) * f_delta(i, delta);
            return a;
        }

        TYPE
        ddelta(TYPE delta, TYPE tau) const
        {
            TYPE da_ddelta = 0.;
            for (std::size_t i = 0; i < this->n.size(); i++) {
                da_ddelta += this->n[i] * f_tau(i, tau) * delta_df_delta_ddelta(i, delta) / delta;
            }
            return da_ddelta;
        }

        TYPE
        dtau(TYPE delta, TYPE tau) const
        {
            TYPE da_dtau = 0.;
            for (std::size_t i = 0; i < this->n.size(); i++) {
                da_dtau += this->n[i] * f_delta(i, delta) * tau_df_tau_dtau(i, tau) / tau;
            }
            return da_dtau;
        }

        TYPE
        d2delta(TYPE delta, TYPE tau) const
        {
            TYPE d2a_ddelta2 = 0.;
            for (std::size_t i = 0; i < this->n.size(); i++)
                d2a_ddelta2 += this->n[i] * f_tau(i, tau) * delta2_d2f_delta_ddelta2(i, delta) /
                               math::pow<2>(delta);
            return d2a_ddelta2;
        }

        TYPE
        d2tau(TYPE delta, TYPE tau) const
        {
            TYPE d2a_dtau2 = 0.;
            for (std::size_t i = 0; i < this->n.size(); i++)
                d2a_dtau2 +=
                    this->n[i] * f_delta(i, delta) * tau2_d2f_tau_dtau2(i, tau) / math::pow<2>(tau);
            return d2a_dtau2;
        }

        TYPE
        d2deltatau(TYPE delta, TYPE tau) const
        {
            TYPE d2a_ddeltatau = 0.;
            for (std::size_t i = 0; i < this->n.size(); i++)
                d2a_ddeltatau += this->n[i] * tau_df_tau_dtau(i, tau) *
                                 delta_df_delta_ddelta(i, delta) / tau / delta;
            return d2a_ddeltatau;
        }

    private:
        TYPE
        f_tau(std::size_t i, double tau) const
        {
            return math::pow(tau, this->t[i]) *
                   std::exp(1.0 /
                            (this->b[i] + this->beta[i] * math::pow<2>(-this->gamma[i] + tau)));
        }

        TYPE
        f_delta(std::size_t i, double delta) const
        {
            return math::pow(delta, this->d[i]) *
                   std::exp(this->eta[i] * math::pow<2>(delta - this->epsilon[i]));
        }

        TYPE
        delta_df_delta_ddelta(std::size_t i, double delta) const
        {
            return (this->d[i] * math::pow(delta, this->d[i]) +
                    2. * math::pow(delta, this->d[i] + 1) * this->eta[i] *
                        (delta - this->epsilon[i])) *
                   std::exp(this->eta[i] * math::pow<2>(delta - this->epsilon[i]));
        }

        TYPE
        tau_df_tau_dtau(std::size_t i, double tau) const
        {
            return (2. * this->beta[i] * math::pow(tau, this->t[i] + 1) * (this->gamma[i] - tau) +
                    this->t[i] * math::pow(tau, this->t[i]) *
                        math::pow<2>(this->b[i] +
                                     this->beta[i] * math::pow<2>(this->gamma[i] - tau))) *
                   std::exp(1.0 /
                            (this->b[i] + this->beta[i] * math::pow<2>(this->gamma[i] - tau))) /
                   math::pow<2>(this->b[i] + this->beta[i] * math::pow<2>(this->gamma[i] - tau));
        }

        TYPE
        delta2_d2f_delta_ddelta2(std::size_t i, double delta) const
        {
            return math::pow(delta, this->d[i]) *
                   (4. * this->d[i] * delta * this->eta[i] * (delta - this->epsilon[i]) +
                    this->d[i] * (this->d[i] - 1) +
                    2. * math::pow<2>(delta) * this->eta[i] *
                        (2. * this->eta[i] * math::pow<2>(delta - this->epsilon[i]) + 1)) *
                   std::exp(this->eta[i] * math::pow<2>(delta - this->epsilon[i]));
        }

        TYPE
        tau2_d2f_tau_dtau2(std::size_t i, double tau) const
        {
            return math::pow(tau, this->t[i]) *
                   (4 * this->beta[i] * this->t[i] * tau *
                        math::pow<2>(this->b[i] +
                                     this->beta[i] * math::pow<2>(this->gamma[i] - tau)) *
                        (this->gamma[i] - tau) +
                    2 * this->beta[i] * math::pow<2>(tau) *
                        (4 * this->beta[i] *
                             (this->b[i] + this->beta[i] * math::pow<2>(this->gamma[i] - tau)) *
                             math::pow<2>(this->gamma[i] - tau) +
                         2 * this->beta[i] * math::pow<2>(this->gamma[i] - tau) -
                         math::pow<2>(this->b[i] +
                                      this->beta[i] * math::pow<2>(this->gamma[i] - tau))) +
                    this->t[i] *
                        math::pow(this->b[i] + this->beta[i] * math::pow<2>(this->gamma[i] - tau),
                                  4) *
                        (this->t[i] - 1)) *
                   std::exp(1.0 /
                            (this->b[i] + this->beta[i] * math::pow<2>(this->gamma[i] - tau))) /
                   math::pow(this->b[i] + this->beta[i] * math::pow<2>(this->gamma[i] - tau), 4);
        }

        std::vector<double> n;
        std::vector<double> d;
        std::vector<double> t;
        std::vector<double> eta;
        std::vector<double> beta;
        std::vector<double> gamma;
        std::vector<double> epsilon;
        std::vector<double> b;
    };

protected:
    double delta(double rho) const;
    double tau(double T) const;

    double temperature(double u, double tau, double da_dt) const;
    double pressure(double rho, double T, double delta, double da_dd) const;
    double internal_energy(double T, double tau, double da_dt) const;
    double enthalphy(double T, double delta, double tau, double da_dt, double da_dd) const;
    double sound_speed(double T,
                       double delta,
                       double tau,
                       double da_dd,
                       double d2a_dd2,
                       double d2a_ddt,
                       double d2a_dt2) const;
    double entropy(double tau, double a, double da_dt) const;
    double heat_capacity_isobaric(double delta,
                                  double tau,
                                  double da_dd,
                                  double d2a_dt2,
                                  double d2a_dd2,
                                  double d2a_ddt) const;
    double heat_capacity_isochoric(double tau, double d2a_dt2) const;
};

} // namespace fprops
