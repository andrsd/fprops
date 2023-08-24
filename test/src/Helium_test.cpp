#include "gtest/gtest.h"
#include "Helium.h"

using namespace fprops;

namespace {

// T = 280 K, p = 1 MPa
SinglePhaseFluidProperties::Props gold1 = { 877864.48974420107,
                                            0.584485817263523,
                                            1.7109055009783694,
                                            1.e6,
                                            280.0,
                                            1.9050487777800348e-05,
                                            5193.7491442602677,
                                            3118.4171250330392,
                                            22983.561984600798,
                                            0.14939414616319868,
                                            1462350.3070077244,
                                            989.0570220636348 };

} // namespace

TEST(HeliumTest, rho_T)
{
    Helium fp;

    double rho = 1.7109055009783694;
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
    EXPECT_DOUBLE_EQ(props.v, gold1.v);
    EXPECT_DOUBLE_EQ(props.s, gold1.s);
    EXPECT_DOUBLE_EQ(props.h, gold1.h);
    EXPECT_DOUBLE_EQ(props.w, gold1.w);
}

TEST(HeliumTest, rho_p)
{
    Helium fp;

    double rho = 1.7109055009783694;
    double p = 1.e6;
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

TEST(HeliumTest, p_T)
{
    Helium fp;

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
    EXPECT_DOUBLE_EQ(props.v, gold1.v);
    EXPECT_DOUBLE_EQ(props.s, gold1.s);
    EXPECT_DOUBLE_EQ(props.h, gold1.h);
    EXPECT_DOUBLE_EQ(props.w, gold1.w);
}

TEST(HeliumTest, v_u)
{
    Helium fp;

    double T = 280.0;
    double p = 1.0e6;
    SinglePhaseFluidProperties::Props state0 = fp.p_T(p, T);

    EXPECT_THROW(auto f = fp.v_u(state0.v, state0.u), std::runtime_error);
    /*
        SinglePhaseFluidProperties::Props props = fp.v_u(state0.v, state0.u);

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
    */
}
