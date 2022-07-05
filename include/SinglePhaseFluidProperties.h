#pragma once

#include "FluidProperties.h"

namespace fprops {

/// Base class for single phase fluids
///
class SinglePhaseFluidProperties : public FluidProperties {
public:
    SinglePhaseFluidProperties();
};

} // namespace fprops
