#include "gtest/gtest.h"
#include "Helium.h"

using namespace fprops;

TEST(HeliumTest, p_T)
{
    Helium fp;

    // 280 K, 1 MPa
    double T = 280.0;
    double p = 1.0e6;
    SinglePhaseFluidProperties::Props props = fp.p_T(p, T);

    EXPECT_DOUBLE_EQ(props.p, p);
    EXPECT_DOUBLE_EQ(props.T, T);

    EXPECT_DOUBLE_EQ(props.rho, 1.7109055009783694);
    EXPECT_DOUBLE_EQ(props.h, 1462350.3070077244);
    EXPECT_DOUBLE_EQ(props.u, 877864.48974420107);

    EXPECT_DOUBLE_EQ(props.s, 22983.561984600798);
    EXPECT_DOUBLE_EQ(props.cp, 5193.7491442602677);
    EXPECT_DOUBLE_EQ(props.cv, 3118.4171250330392);
    EXPECT_DOUBLE_EQ(props.w, 989.0570220636348);

    EXPECT_DOUBLE_EQ(props.mu, 1.9050487777800348e-05);
    EXPECT_DOUBLE_EQ(props.k, 0.14939414616319868);
}

TEST(HeliumTest, v_u)
{
    Helium fp;

    // 280 K, 1 MPa
    double T = 280.0;
    double p = 1.0e6;
    SinglePhaseFluidProperties::Props state0 = fp.p_T(p, T);

    EXPECT_THROW(auto f = fp.v_u(state0.v, state0.u), std::runtime_error);
    /*
        SinglePhaseFluidProperties::Props props = fp.v_u(state0.v, state0.u);

        EXPECT_DOUBLE_EQ(props.p, p);
        EXPECT_DOUBLE_EQ(props.T, T);
        EXPECT_DOUBLE_EQ(props.rho, 1.7109055009783694);
        EXPECT_DOUBLE_EQ(props.h, 1462350.3070077244);
        EXPECT_DOUBLE_EQ(props.u, 877864.48974420107);
        EXPECT_DOUBLE_EQ(props.s, 22983.561984600798);
        EXPECT_DOUBLE_EQ(props.cp, 5193.7491442602677);
        EXPECT_DOUBLE_EQ(props.cv, 3118.4171250330392);
        EXPECT_DOUBLE_EQ(props.w, 989.0570220636348);

        EXPECT_DOUBLE_EQ(props.mu, 1.7090507109297636e-05);
        EXPECT_DOUBLE_EQ(props.k, 0.024857419011067187);
    */
}
