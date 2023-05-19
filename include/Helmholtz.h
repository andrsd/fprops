#pragma once

#include "SinglePhaseFluidProperties.h"
#include "Numerics.h"
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

    [[nodiscard]] Props p_T(double p, double T) const override;
    [[nodiscard]] Props v_u(double v, double u) const override;

protected:
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
                sum += this->n[i] * std::pow(tau, this->t[i]);
            return sum;
        }

        T
        dtau(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < n.size(); ++i)
                sum += this->n[i] * this->t[i] * std::pow(tau, this->t[i] - 1);
            return sum;
        }

        T
        d2tau(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); ++i)
                sum += this->n[i] * this->t[i] * (this->t[i] - 1) * std::pow(tau, this->t[i] - 2);
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
                sum += this->n[i] * sqr(this->theta[i]) * this->c[i] * this->d[i] * exp_theta_tau /
                       sqr(this->c[i] + this->d[i] * exp_theta_tau);
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
                sum += this->n[i] * std::pow(delta, this->d[i]) * std::pow(tau, this->t[i]);
            return sum;
        }

        T
        ddelta(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum += this->n[i] * this->d[i] * std::pow(delta, this->d[i] - 1) *
                       std::pow(tau, this->t[i]);
            return sum;
        }

        T
        dtau(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum += this->n[i] * std::pow(delta, this->d[i]) * this->t[i] *
                       std::pow(tau, this->t[i] - 1);
            return sum;
        }

        T
        d2delta(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum += this->n[i] * this->d[i] * (this->d[i] - 1) *
                       std::pow(delta, this->d[i] - 2) * std::pow(tau, this->t[i]);
            return sum;
        }

        T
        d2tau(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum += this->n[i] * std::pow(delta, this->d[i]) * this->t[i] * (this->t[i] - 1) *
                       std::pow(tau, this->t[i] - 2);
            return sum;
        }

        T
        d2deltatau(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum += this->n[i] * this->d[i] * std::pow(delta, this->d[i] - 1) * this->t[i] *
                       std::pow(tau, this->t[i] - 1);
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
                    (sqr(this->l[i]) * std::pow(delta, 2 * this->l[i]) +
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
                       std::exp(-this->eta[i] * sqr(delta - this->epsilon[i]) -
                                this->beta[i] * sqr(tau - this->gamma[i]));
            return sum;
        }

        T
        ddelta(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum += this->n[i] * std::pow(delta, this->d[i] - 1) * std::pow(tau, this->t[i]) *
                       (this->d[i] - 2.0 * delta * this->eta[i] * (delta - this->epsilon[i])) *
                       std::exp(-this->eta[i] * sqr(delta - this->epsilon[i]) -
                                this->beta[i] * sqr(tau - this->gamma[i]));
            return sum;
        }

        T
        dtau(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum += this->n[i] * std::pow(delta, this->d[i]) * std::pow(tau, this->t[i] - 1) *
                       (2.0 * this->beta[i] * (this->gamma[i] - tau) * tau + this->t[i]) *
                       std::exp(-this->eta[i] * sqr(delta - this->epsilon[i]) -
                                this->beta[i] * sqr(tau - this->gamma[i]));
            return sum;
        }

        T
        d2delta(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum += this->n[i] * std::pow(delta, this->d[i] - 2) * std::pow(tau, this->t[i]) *
                       (2 * sqr(delta) * this->eta[i] *
                            (2 * this->eta[i] * sqr(delta - this->epsilon[i]) - 1) +
                        sqr(this->d[i]) +
                        this->d[i] * (4 * delta * this->eta[i] * (this->epsilon[i] - delta) - 1)) *
                       std::exp(-this->eta[i] * sqr(delta - this->epsilon[i]) -
                                this->beta[i] * sqr(tau - this->gamma[i]));
            return sum;
        }

        T
        d2tau(T delta, T tau) const
        {
            T sum = 0;
            for (std::size_t i = 0; i < this->n.size(); i++)
                sum += this->n[i] * std::pow(delta, this->d[i]) * std::pow(tau, this->t[i] - 2) *
                       (2 * this->beta[i] * sqr(tau) *
                            (2 * this->beta[i] * sqr(this->gamma[i] - tau) - 1) +
                        sqr(this->t[i]) +
                        this->t[i] * (4 * this->beta[i] * (this->gamma[i] - tau) * tau - 1)) *
                       std::exp(-this->eta[i] * sqr(delta - this->epsilon[i]) -
                                this->beta[i] * sqr(tau - this->gamma[i]));
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
                       std::exp(-this->beta[i] * sqr(this->gamma[i] - tau) -
                                this->eta[i] * sqr(delta - this->epsilon[i]));
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
};

} // namespace fprops
