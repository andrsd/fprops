#include "gtest/gtest.h"
#include "CarbonDioxide.h"

using namespace fprops;

namespace {

// T = 280 K, p = 1 MPa
SinglePhaseFluidProperties::Props gold1 = { 430888.05580697849,
                                            0.049506644012416195,
                                            20.199309000812121,
                                            1e6,
                                            280,
                                            1.4150532169060509e-05,
                                            925.17988930107481,
                                            670.91985675762839,
                                            2225.7418437948568,
                                            0.015727797767537487,
                                            480394.69981939468,
                                            252.32654137014907 };

} // namespace

TEST(CarbonDioxide, rho_T)
{
    CarbonDioxide fp;

    double rho = 20.199309000812121;
    double T = 280;
    SinglePhaseFluidProperties::Props props = fp.rho_T(rho, T);

    EXPECT_DOUBLE_EQ(props.rho, gold1.rho);
    EXPECT_DOUBLE_EQ(props.T, gold1.T);
    EXPECT_DOUBLE_EQ(props.p, gold1.p);
    EXPECT_DOUBLE_EQ(props.u, gold1.u);
    EXPECT_DOUBLE_EQ(props.cv, gold1.cv);
    EXPECT_DOUBLE_EQ(props.cp, gold1.cp);
    EXPECT_DOUBLE_EQ(props.mu, gold1.mu);
    EXPECT_DOUBLE_EQ(props.k, gold1.k);
    EXPECT_DOUBLE_EQ(props.v, gold1.v);
    EXPECT_DOUBLE_EQ(props.s, gold1.s);
    EXPECT_DOUBLE_EQ(props.h, gold1.h);
    EXPECT_DOUBLE_EQ(props.w, gold1.w);
}

TEST(CarbonDioxide, rho_p)
{
    CarbonDioxide fp;

    double rho = 20.199309000812121;
    double p = 1e6;
    SinglePhaseFluidProperties::Props props = fp.rho_p(rho, p);

    EXPECT_DOUBLE_EQ(props.rho, gold1.rho);
    EXPECT_DOUBLE_EQ(props.T, gold1.T);
    EXPECT_DOUBLE_EQ(props.p, gold1.p);
    EXPECT_DOUBLE_EQ(props.u, gold1.u);
    EXPECT_DOUBLE_EQ(props.cv, gold1.cv);
    EXPECT_DOUBLE_EQ(props.cp, gold1.cp);
    EXPECT_DOUBLE_EQ(props.mu, gold1.mu);
    EXPECT_DOUBLE_EQ(props.k, gold1.k);
    EXPECT_DOUBLE_EQ(props.v, gold1.v);
    EXPECT_DOUBLE_EQ(props.s, gold1.s);
    EXPECT_DOUBLE_EQ(props.h, gold1.h);
    EXPECT_DOUBLE_EQ(props.w, gold1.w);
}

TEST(CarbonDioxide, p_T)
{
    CarbonDioxide fp;

    double T = 280;
    double p = 1e6;
    SinglePhaseFluidProperties::Props props = fp.p_T(p, T);

    EXPECT_DOUBLE_EQ(props.rho, gold1.rho);
    EXPECT_DOUBLE_EQ(props.T, gold1.T);
    EXPECT_DOUBLE_EQ(props.p, gold1.p);
    EXPECT_DOUBLE_EQ(props.u, gold1.u);
    EXPECT_DOUBLE_EQ(props.cv, gold1.cv);
    EXPECT_DOUBLE_EQ(props.cp, gold1.cp);
    EXPECT_DOUBLE_EQ(props.mu, gold1.mu);
    EXPECT_DOUBLE_EQ(props.k, gold1.k);
    EXPECT_DOUBLE_EQ(props.v, gold1.v);
    EXPECT_DOUBLE_EQ(props.s, gold1.s);
    EXPECT_DOUBLE_EQ(props.h, gold1.h);
    EXPECT_DOUBLE_EQ(props.w, gold1.w);
}

TEST(CarbonDioxide, v_u)
{
    CarbonDioxide fp;

    double v = 0.049506644012416195;
    double u = 430888.05580697849;
    SinglePhaseFluidProperties::Props props = fp.v_u(v, u);

    EXPECT_DOUBLE_EQ(props.rho, gold1.rho);
    EXPECT_DOUBLE_EQ(props.T, gold1.T);
    EXPECT_DOUBLE_EQ(props.p, gold1.p);
    EXPECT_DOUBLE_EQ(props.u, gold1.u);
    EXPECT_DOUBLE_EQ(props.cv, gold1.cv);
    EXPECT_DOUBLE_EQ(props.cp, gold1.cp);
    EXPECT_DOUBLE_EQ(props.mu, gold1.mu);
    EXPECT_DOUBLE_EQ(props.k, gold1.k);
    EXPECT_DOUBLE_EQ(props.v, gold1.v);
    EXPECT_DOUBLE_EQ(props.s, gold1.s);
    EXPECT_DOUBLE_EQ(props.h, gold1.h);
    EXPECT_DOUBLE_EQ(props.w, gold1.w);
}
