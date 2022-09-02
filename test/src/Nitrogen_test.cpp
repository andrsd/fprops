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

    EXPECT_DOUBLE_EQ(props.rho, 12.074993450711256);
    EXPECT_DOUBLE_EQ(props.h, 288083.48983915086);
    EXPECT_DOUBLE_EQ(props.u, 205267.70993153556);
    EXPECT_DOUBLE_EQ(props.s, 6083.1854964551458);
    EXPECT_DOUBLE_EQ(props.cp, 1058.6154679415292);
    EXPECT_DOUBLE_EQ(props.cv, 745.56949807601234);
    EXPECT_DOUBLE_EQ(props.w, 342.35848440003775);

    EXPECT_DOUBLE_EQ(props.mu, 1.709050710929316e-05);
    EXPECT_DOUBLE_EQ(props.k, 0.024857419011055305);
}

TEST(NitrogenTest, v_u)
{
    Nitrogen fp;

    // 280 K, 1 MPa
    double T = 280.0;
    double p = 1.0e6;
    SinglePhaseFluidProperties::Props state0 = fp.p_T(p, T);

    SinglePhaseFluidProperties::Props props = fp.v_u(state0.v, state0.u);

    EXPECT_DOUBLE_EQ(props.p, p);
    EXPECT_DOUBLE_EQ(props.T, T);
    EXPECT_DOUBLE_EQ(props.rho, 12.074993450711256);
    EXPECT_DOUBLE_EQ(props.h, 288083.48983915086);
    EXPECT_DOUBLE_EQ(props.u, 205267.70993153556);
    EXPECT_DOUBLE_EQ(props.s, 6083.1854964551458);
    EXPECT_DOUBLE_EQ(props.cp, 1058.6154679415292);
    EXPECT_DOUBLE_EQ(props.cv, 745.56949807601234);
    EXPECT_DOUBLE_EQ(props.w, 342.35848440003775);

    EXPECT_DOUBLE_EQ(props.mu, 1.709050710929316e-05);
    EXPECT_DOUBLE_EQ(props.k, 0.024857419011055305);
}
