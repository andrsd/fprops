#include "gmock/gmock.h"
#include "fprops/fprops.h"
#include "fprops/Air.h"
#include "fprops/CarbonDioxide.h"
#include "fprops/Helium.h"
#include "fprops/Nitrogen.h"

using namespace fprops;
using namespace testing;

TEST(FpropsTest, from_name)
{
    auto air = fprops::from_name("air");
    EXPECT_THAT(dynamic_cast<Air *>(air), NotNull());

    auto co2 = fprops::from_name("co2");
    EXPECT_THAT(dynamic_cast<CarbonDioxide *>(co2), NotNull());

    auto he = fprops::from_name("helium");
    EXPECT_THAT(dynamic_cast<Helium *>(he), NotNull());

    auto n2 = fprops::from_name("nitrogen");
    EXPECT_THAT(dynamic_cast<Nitrogen *>(n2), NotNull());
}
