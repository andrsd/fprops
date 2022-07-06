#pragma once

#include "FluidProperties.h"

namespace fprops {

/// Base class for single phase fluids
///
class SinglePhaseFluidProperties : public FluidProperties {
public:
    SinglePhaseFluidProperties();

    /// Computed properties referring to a thermodynamical state
    struct Props {
        /// Internal energy
        double u;
        /// Specific volume
        double v;
        /// Density
        double rho;
        /// Pressure
        double p;
        /// Temperature
        double T;
        /// Dynamic viscosity
        double mu;
        /// Isobaric heat capacity
        double cp;
        /// Isochoric heat capacity
        double cv;
        /// Entropy
        double s;
        /// Thermal conductivity
        double k;
        /// Specific enthalpy
        double h;
        /// Speed of sound
        double w;

        Props();
    };

    /// Compute thermodynamical state given pressure and temperature
    ///
    /// @param p Pressure
    /// @param T Temperature
    virtual Props p_T(double p, double T) = 0;

    /// Compute thermodynamical state given specific volume and internal energy
    ///
    /// @param v Specific volume
    /// @param u Internal energy
    virtual Props v_u(double v, double u) = 0;
};

} // namespace fprops
