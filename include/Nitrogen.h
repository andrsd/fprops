#pragma once

#include "Helmholtz.h"

namespace fprops {

/// Nitrogen (N2) fluid properties
///
/// Thermodynamic properties calculated from:
/// Span,. Lemmon, Jacobsen, Wagner and Yokozeki, A reference equation of state for the
/// thermodynamic properties of nitrogen for temperatures from 63.151 to 1000 K and pressures to
/// 2200 MPa, Journal of Physical and Chemical Reference Data, 29, 1361--1433 (2000)
class Nitrogen : public Helmholtz {
public:
    Nitrogen();

protected:
    virtual double alpha(double delta, double tau) override;
    virtual double dalpha_ddelta(double delta, double tau) override;
    virtual double dalpha_dtau(double delta, double tau) override;
    virtual double d2alpha_ddelta2(double delta, double tau) override;
    virtual double d2alpha_dtau2(double delta, double tau) override;
    virtual double d2alpha_ddeltatau(double delta, double tau) override;
};

} // namespace fprops
