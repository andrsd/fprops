#pragma once

#include "FluidProperties.h"

namespace fprops {

/// Base class for single phase fluids
///
class SinglePhaseFluidProperties : public FluidProperties {
public:
    SinglePhaseFluidProperties();
    virtual ~SinglePhaseFluidProperties();

    /// Computed properties referring to a thermodynamical state
    struct Props {
        /// Internal energy \f$[J/kg]\f$
        double u;
        /// Specific volume \f$[m^3/kg]\f$
        double v;
        /// Density \f$[kg/m^3]\f$
        double rho;
        /// Pressure \f$[Pa]\f$
        double p;
        /// Temperature \f$[K]\f$
        double T;
        /// Dynamic viscosity \f$[Pa-s]\f$
        double mu;
        /// Isobaric heat capacity \f$[J/(kg-K)]\f$
        double cp;
        /// Isochoric heat capacity \f$[J/(kg-K)]\f$
        double cv;
        /// Entropy \f$[J/(kg-K)]\f$
        double s;
        /// Thermal conductivity \f$[W/(m-K)]\f$
        double k;
        /// Specific enthalpy \f$[J/kg]\f$
        double h;
        /// Speed of sound \f$[m/s]\f$
        double w;

        Props();
    };

    /// Compute thermodynamical state given pressure and temperature
    ///
    /// @param p Pressure \f$[Pa]\f$
    /// @param T Temperature \f$[K]\f$
    virtual Props p_T(double p, double T) const = 0;

    /// Compute thermodynamical state given specific volume and internal energy
    ///
    /// @param v Specific volume \f$[m^3/kg]\f$
    /// @param u Specific internal energy \f$[J/kg]\f$
    virtual Props v_u(double v, double u) const = 0;
};

} // namespace fprops
