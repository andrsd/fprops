#include "gtest/gtest.h"
#include "CarbonDioxide.h"

using namespace fprops;

TEST(CarbonDioxide, p_T)
{
    CarbonDioxide fp;

    double T = 280;
    double p = 1e6;
    SinglePhaseFluidProperties::Props props = fp.p_T(p, T);

    EXPECT_DOUBLE_EQ(props.p, p);
    EXPECT_DOUBLE_EQ(props.T, T);
    EXPECT_DOUBLE_EQ(props.rho, 20.199309000812121);
    EXPECT_DOUBLE_EQ(props.u, 430888.05580697849);
    EXPECT_DOUBLE_EQ(props.cv, 670.91985675762839);
    EXPECT_DOUBLE_EQ(props.cp, 925.17988930107481);
    EXPECT_DOUBLE_EQ(props.v, 0.049506644012416195);
    EXPECT_DOUBLE_EQ(props.s, 2225.7418437948568);
    EXPECT_DOUBLE_EQ(props.h, 480394.69981939468);
    EXPECT_DOUBLE_EQ(props.w, 252.32654137014907);
    EXPECT_DOUBLE_EQ(props.mu, 1.4150532169060509e-05);
    EXPECT_DOUBLE_EQ(props.k, 0.015727797767537487);
}

TEST(CarbonDioxide, v_u)
{
    CarbonDioxide fp;

    double v = 0.049506644012416195;
    double u = 430888.05580697849;
    SinglePhaseFluidProperties::Props props = fp.v_u(v, u);

    EXPECT_DOUBLE_EQ(props.v, v);
    EXPECT_DOUBLE_EQ(props.u, u);
    EXPECT_DOUBLE_EQ(props.p, 1e6);
    EXPECT_DOUBLE_EQ(props.T, 280);
    EXPECT_DOUBLE_EQ(props.rho, 20.199309000812121);
    EXPECT_DOUBLE_EQ(props.cv, 670.91985675762839);
    EXPECT_DOUBLE_EQ(props.cp, 925.17988930107481);
    EXPECT_DOUBLE_EQ(props.s, 2225.7418437948568);
    EXPECT_DOUBLE_EQ(props.h, 480394.69981939468);
    EXPECT_DOUBLE_EQ(props.w, 252.32654137014907);
    EXPECT_DOUBLE_EQ(props.mu, 1.4150532169060509e-05);
    EXPECT_DOUBLE_EQ(props.k, 0.015727797767537487);
}
