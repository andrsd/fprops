#include "gtest/gtest.h"
#include "Nitrogen.h"

using namespace fprops;

TEST(NitrogenTest, p_T)
{
    Nitrogen fp;

    // 280 K, 1 MPa
    double T = 280.0;
    double p = 1.0e6;
    SinglePhaseFluidProperties::Props props = fp.p_T(p, T);

    EXPECT_DOUBLE_EQ(props.p, p);
    EXPECT_DOUBLE_EQ(props.T, T);

    EXPECT_DOUBLE_EQ(props.rho, 12.074993445601923);
    EXPECT_DOUBLE_EQ(props.h, 288083.49016938655);
    EXPECT_DOUBLE_EQ(props.u, 205267.7102267291);
    EXPECT_DOUBLE_EQ(props.s, 6083.1855363670165);
    EXPECT_DOUBLE_EQ(props.cp, 1058.6154654154957);
    EXPECT_DOUBLE_EQ(props.cv, 745.56949769606433);
    EXPECT_DOUBLE_EQ(props.w, 342.35848421182038);
}
