// SPDX-FileCopyrightText: 2022 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#pragma once

#include "fprops/state.h"
#include <memory>

namespace fprops {

/// Single phase fluids
///
class SinglePhaseFluidProperties {
private:
    struct AbstractProperties {
        virtual ~AbstractProperties() = default;
        virtual State rho_T(double rho, double T) const = 0;
        virtual State rho_p(double rho, double p) const = 0;
        virtual State p_T(double p, double T) const = 0;
        virtual State v_u(double v, double u) const = 0;
        virtual State h_s(double h, double s) const = 0;
    };

    template <typename FP>
    struct Properties : public AbstractProperties {
        FP fprops_;

        explicit Properties(FP fprops) : fprops_(fprops) {}

        State
        rho_T(double rho, double T) const override
        {
            return this->fprops_.rho_T(rho, T);
        }

        State
        rho_p(double rho, double p) const override
        {
            return this->fprops_.rho_p(rho, p);
        }

        State
        p_T(double p, double T) const override
        {
            return this->fprops_.p_T(p, T);
        }

        State
        v_u(double v, double u) const override
        {
            return this->fprops_.v_u(v, u);
        }

        State
        h_s(double h, double s) const override
        {
            return this->fprops_.h_s(h, s);
        }
    };

    std::shared_ptr<AbstractProperties> impl_;

public:
    SinglePhaseFluidProperties() = default;

    // constructor from raw pointer to an implementation
    template <typename FPROPS>
    explicit SinglePhaseFluidProperties(FPROPS fprops) :
        impl_(std::make_shared<Properties<FPROPS>>(fprops))
    {
    }

    operator bool() const { return this->impl_ != nullptr; }

    /// Compute thermodynamical state given density and temperature
    ///
    /// @param rho Density \f$[kg/m^3]\f$
    /// @param T Temperature \f$[K]\f$
    [[nodiscard]] State rho_T(double rho, double T) const;

    /// Compute thermodynamical state given density and pressure
    ///
    /// @param rho Density \f$[kg/m^3]\f$
    /// @param p Pressure \f$[Pa]\f$
    [[nodiscard]] State rho_p(double rho, double p) const;

    /// Compute thermodynamical state given pressure and temperature
    ///
    /// @param p Pressure \f$[Pa]\f$
    /// @param T Temperature \f$[K]\f$
    [[nodiscard]] State p_T(double p, double T) const;

    /// Compute thermodynamical state given specific volume and internal energy
    ///
    /// @param v Specific volume \f$[m^3/kg]\f$
    /// @param u Specific internal energy \f$[J/kg]\f$
    [[nodiscard]] State v_u(double v, double u) const;

    /// Compute thermodynamical state given specific enthalpy and entropy
    ///
    /// @param h Specific enthalpy \f$[J/kg]\f$
    /// @param s Entropy \f$[J/(kg-K)]\f$
    [[nodiscard]] State h_s(double h, double s) const;

public:
    /// Construct fluid properties from a fluid name
    ///
    /// @param name Fluid name
    /// @return Fluid properties
    static SinglePhaseFluidProperties from_name(const std::string & name);
};

} // namespace fprops
