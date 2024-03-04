#pragma once

#include <string>

namespace fprops {

/// Computed properties referring to a thermodynamical state
struct State {
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

    State();
    /// Constructor
    ///
    /// @param u Internal energy \f$[J/kg]\f$
    /// @param v Specific volume \f$[m^3/kg]\f$
    /// @param rho Density \f$[kg/m^3]\f$
    /// @param p Pressure \f$[Pa]\f$
    /// @param T Temperature \f$[K]\f$
    /// @param mu Dynamic viscosity \f$[Pa-s]\f$
    /// @param cp Isobaric heat capacity \f$[J/(kg-K)]\f$
    /// @param cv Isochoric heat capacity \f$[J/(kg-K)]\f$
    /// @param s Entropy \f$[J/(kg-K)]\f$
    /// @param k Thermal conductivity \f$[W/(m-K)]\f$
    /// @param h Specific enthalpy \f$[J/kg]\f$
    /// @param w Speed of sound \f$[m/s]\f$
    State(double u,
          double v,
          double rho,
          double p,
          double T,
          double mu,
          double cp,
          double cv,
          double s,
          double k,
          double h,
          double w);

    [[nodiscard]] auto to_string() const -> std::string;
};

} // namespace fprops
