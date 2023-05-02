#include "Air.h"
#include "Numerics.h"
#include <cmath>

namespace fprops {

Air::Air() : Helmholtz(8.314510, 28.9586e-3, 342.60340488, 132.5306) {}

double
Air::alpha(double delta, double tau) const
{
    double alpha0 = 0;
    double alphar = 0.0;
    return alpha0 + alphar;
}

double
Air::dalpha_ddelta(double delta, double tau) const
{
    double dalpha0 = 0;
    double dalphar = 0.0;
    return dalpha0 + dalphar;
}

double
Air::dalpha_dtau(double delta, double tau) const
{
    const double dalpha0 = 0;
    double dalphar = 0.0;
    return dalpha0 + dalphar;
}

double
Air::d2alpha_ddelta2(double delta, double tau) const
{
    const double dalpha0 = 0;
    double dalphar = 0.0;
    return dalpha0 + dalphar;
}

double
Air::d2alpha_dtau2(double delta, double tau) const
{
    double dalpha0 = 0.;
    double dalphar = 0.0;
    return dalpha0 + dalphar;
}

double
Air::d2alpha_ddeltatau(double delta, double tau) const
{
    double dalphar = 0.0;
    return dalphar;
}

double
Air::mu_from_rho_T(double rho, double T) const
{
    double mu0 = 0.0;
    double mur = 0.0;
    return (mu0 + mur) * 1.0e-6;
}

double
Air::k_from_rho_T(double rho, double T) const
{
    double lambda0 = 0.0;
    double lambdar = 0.0;
    return (lambda0 + lambdar) * 1.0e-3;
}

} // namespace fprops
