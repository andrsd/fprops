#include "gmock/gmock.h"
#include "fprops/single_phase_fluid_properties.h"

using namespace fprops;
using namespace testing;

TEST(FpropsTest, from_name)
{
    auto air = SinglePhaseFluidProperties::from_name("air");
    // TODO: check that it is air

    auto co2 = SinglePhaseFluidProperties::from_name("co2");
    // TODO: check that it is co2

    auto he = SinglePhaseFluidProperties::from_name("helium");
    // TODO: check that it is helium

    auto n2 = SinglePhaseFluidProperties::from_name("nitrogen");
    // TODO: check that it is nitrogen
}
