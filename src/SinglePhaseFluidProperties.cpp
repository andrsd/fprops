#include "SinglePhaseFluidProperties.h"

namespace fprops {

SinglePhaseFluidProperties::Props::Props()
{
    this->u = 0.;
    this->v = 0.;
    this->rho = 0.;
    this->p = 0.;
    this->T = 0.;
    this->mu = 0.;
    this->cp = 0.;
    this->cv = 0.;
    this->s = 0.;
    this->k = 0.;
    this->h = 0.;
    this->w = 0.;
}

SinglePhaseFluidProperties::SinglePhaseFluidProperties() : FluidProperties() {}

} // namespace fprops
