#include "gtest/gtest.h"
#include "Nitrogen.h"

using namespace fprops;

namespace {

// T = 280 K, p = 1 MPa
SinglePhaseFluidProperties::Props gold1 = { 205267.70993394594,
                                            0.082815779905282,
                                            12.074993451051515,
                                            1.0e6,
                                            280.0,
                                            1.7090507109297636e-05,
                                            1058.6154681901673,
                                            745.56949823705611,
                                            6083.1854964583363,
                                            0.024857419011067187,
                                            288083.48983922758,
                                            342.35848437431741 };

} // namespace

TEST(NitrogenTest, rho_T)
{
    Nitrogen fp;

    double rho = 12.074993451051515;
    double T = 280.0;
    SinglePhaseFluidProperties::Props props = fp.rho_T(rho, T);

    EXPECT_DOUBLE_EQ(props.rho, gold1.rho);
    EXPECT_DOUBLE_EQ(props.T, gold1.T);
    EXPECT_DOUBLE_EQ(props.p, gold1.p);
    EXPECT_DOUBLE_EQ(props.u, gold1.u);
    EXPECT_DOUBLE_EQ(props.cv, gold1.cv);
    EXPECT_DOUBLE_EQ(props.cp, gold1.cp);
    EXPECT_DOUBLE_EQ(props.mu, gold1.mu);
    EXPECT_DOUBLE_EQ(props.k, gold1.k);
    EXPECT_NEAR(props.v, gold1.v, 1e-9);
    EXPECT_DOUBLE_EQ(props.s, gold1.s);
    EXPECT_DOUBLE_EQ(props.h, gold1.h);
    EXPECT_DOUBLE_EQ(props.w, gold1.w);
}

TEST(NitrogenTest, rho_p)
{
    Nitrogen fp;

    double rho = 12.074993451051515;
    double p = 1.0e6;
    SinglePhaseFluidProperties::Props props = fp.rho_p(rho, p);

    EXPECT_DOUBLE_EQ(props.rho, gold1.rho);
    EXPECT_DOUBLE_EQ(props.T, gold1.T);
    EXPECT_DOUBLE_EQ(props.p, gold1.p);
    EXPECT_DOUBLE_EQ(props.u, gold1.u);
    EXPECT_DOUBLE_EQ(props.cv, gold1.cv);
    EXPECT_DOUBLE_EQ(props.cp, gold1.cp);
    EXPECT_DOUBLE_EQ(props.mu, gold1.mu);
    EXPECT_DOUBLE_EQ(props.k, gold1.k);
    EXPECT_NEAR(props.v, gold1.v, 1e-9);
    EXPECT_DOUBLE_EQ(props.s, gold1.s);
    EXPECT_DOUBLE_EQ(props.h, gold1.h);
    EXPECT_DOUBLE_EQ(props.w, gold1.w);
}

TEST(NitrogenTest, p_T)
{
    Nitrogen fp;

    double T = 280.0;
    double p = 1.0e6;
    SinglePhaseFluidProperties::Props props = fp.p_T(p, T);

    EXPECT_DOUBLE_EQ(props.rho, gold1.rho);
    EXPECT_DOUBLE_EQ(props.T, gold1.T);
    EXPECT_DOUBLE_EQ(props.p, gold1.p);
    EXPECT_DOUBLE_EQ(props.u, gold1.u);
    EXPECT_DOUBLE_EQ(props.cv, gold1.cv);
    EXPECT_DOUBLE_EQ(props.cp, gold1.cp);
    EXPECT_DOUBLE_EQ(props.mu, gold1.mu);
    EXPECT_DOUBLE_EQ(props.k, gold1.k);
    EXPECT_NEAR(props.v, gold1.v, 1e-9);
    EXPECT_DOUBLE_EQ(props.s, gold1.s);
    EXPECT_DOUBLE_EQ(props.h, gold1.h);
    EXPECT_DOUBLE_EQ(props.w, gold1.w);
}

TEST(NitrogenTest, v_u)
{
    Nitrogen fp;

    double v = 0.082815779905282;
    double u = 205267.70993394594;
    SinglePhaseFluidProperties::Props props = fp.v_u(v, u);

    EXPECT_NEAR(props.rho, gold1.rho, 1e-9);
    EXPECT_DOUBLE_EQ(props.T, gold1.T);
    EXPECT_NEAR(props.p, gold1.p, 1e-8);
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
