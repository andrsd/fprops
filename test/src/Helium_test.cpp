#include "gtest/gtest.h"
#include "ExceptionTestMacros.h"
#include "fprops/Helium.h"

using namespace fprops;

namespace {

// T = 280 K, p = 1 MPa
State gold1 = { 877864.48974420107,
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
    auto state = fp.rho_T(rho, T);

    EXPECT_DOUBLE_EQ(state.rho, gold1.rho);
    EXPECT_DOUBLE_EQ(state.T, gold1.T);
    EXPECT_DOUBLE_EQ(state.p, gold1.p);
    EXPECT_DOUBLE_EQ(state.u, gold1.u);
    EXPECT_DOUBLE_EQ(state.cv, gold1.cv);
    EXPECT_DOUBLE_EQ(state.cp, gold1.cp);
    EXPECT_DOUBLE_EQ(state.mu, gold1.mu);
    EXPECT_DOUBLE_EQ(state.k, gold1.k);
    EXPECT_DOUBLE_EQ(state.v, gold1.v);
    EXPECT_DOUBLE_EQ(state.s, gold1.s);
    EXPECT_DOUBLE_EQ(state.h, gold1.h);
    EXPECT_DOUBLE_EQ(state.w, gold1.w);
}

TEST(HeliumTest, rho_p)
{
    Helium fp;

    double rho = 1.7109055009783694;
    double p = 1.e6;
    auto state = fp.rho_p(rho, p);

    EXPECT_DOUBLE_EQ(state.rho, gold1.rho);
    EXPECT_DOUBLE_EQ(state.T, gold1.T);
    EXPECT_DOUBLE_EQ(state.p, gold1.p);
    EXPECT_DOUBLE_EQ(state.u, gold1.u);
    EXPECT_DOUBLE_EQ(state.cv, gold1.cv);
    EXPECT_DOUBLE_EQ(state.cp, gold1.cp);
    EXPECT_DOUBLE_EQ(state.mu, gold1.mu);
    EXPECT_DOUBLE_EQ(state.k, gold1.k);
    EXPECT_DOUBLE_EQ(state.v, gold1.v);
    EXPECT_DOUBLE_EQ(state.s, gold1.s);
    EXPECT_DOUBLE_EQ(state.h, gold1.h);
    EXPECT_DOUBLE_EQ(state.w, gold1.w);
}

TEST(HeliumTest, p_T)
{
    Helium fp;

    double T = 280.0;
    double p = 1.0e6;
    auto state = fp.p_T(p, T);

    EXPECT_DOUBLE_EQ(state.rho, gold1.rho);
    EXPECT_DOUBLE_EQ(state.T, gold1.T);
    EXPECT_DOUBLE_EQ(state.p, gold1.p);
    EXPECT_DOUBLE_EQ(state.u, gold1.u);
    EXPECT_DOUBLE_EQ(state.cv, gold1.cv);
    EXPECT_DOUBLE_EQ(state.cp, gold1.cp);
    EXPECT_DOUBLE_EQ(state.mu, gold1.mu);
    EXPECT_DOUBLE_EQ(state.k, gold1.k);
    EXPECT_DOUBLE_EQ(state.v, gold1.v);
    EXPECT_DOUBLE_EQ(state.s, gold1.s);
    EXPECT_DOUBLE_EQ(state.h, gold1.h);
    EXPECT_DOUBLE_EQ(state.w, gold1.w);
}

TEST(HeliumTest, v_u)
{
    Helium fp;

    double T = 280.0;
    double p = 1.0e6;
    auto state0 = fp.p_T(p, T);

    EXPECT_THROW_MSG(auto f = fp.v_u(state0.v, state0.u), "Newton's method failed to converge");
    /*
        State state = fp.v_u(state0.v, state0.u);

        EXPECT_DOUBLE_EQ(state.rho, gold1.rho);
        EXPECT_DOUBLE_EQ(state.T, gold1.T);
        EXPECT_DOUBLE_EQ(state.p, gold1.p);
        EXPECT_DOUBLE_EQ(state.u, gold1.u);
        EXPECT_DOUBLE_EQ(state.cv, gold1.cv);
        EXPECT_DOUBLE_EQ(state.cp, gold1.cp);
        EXPECT_DOUBLE_EQ(state.mu, gold1.mu);
        EXPECT_DOUBLE_EQ(state.k, gold1.k);
        EXPECT_DOUBLE_EQ(state.v, gold1.v);
        EXPECT_DOUBLE_EQ(state.s, gold1.s);
        EXPECT_DOUBLE_EQ(state.h, gold1.h);
        EXPECT_DOUBLE_EQ(state.w, gold1.w);
    */
}
