#include "gtest/gtest.h"
#include <gmock/gmock-matchers.h>
#include "exception_test_macros.h"
#include "fprops/methane.h"

using namespace fprops;
using namespace testing;

namespace {

// T = 280 K, p = 1 MPa

double w = 431.97873645972703;
double rho = 7.043156624753545;
double h = 859902.14463317941;
double s = 5323.7095834448064;
double T = 280;
double p = 1e6;
double cp = 2257.7903212900533;
double cv = 1680.8602830296586;
double u = 717920.35247972247;
double v = 1. / rho;
double mu = 1.0738075482060213e-05;
double k = 0.03232260677413022;

State gold1 = { u, v, rho, p, T, mu, cp, cv, s, k, h, w };

} // namespace

TEST(MethaneTest, rho_T)
{
    Methane fp;
    auto state = fp.rho_T(rho, T);

    EXPECT_DOUBLE_EQ(state.rho, gold1.rho);
    EXPECT_DOUBLE_EQ(state.T, gold1.T);
    EXPECT_DOUBLE_EQ(state.p, gold1.p);
    EXPECT_DOUBLE_EQ(state.u, gold1.u);
    EXPECT_DOUBLE_EQ(state.cv, gold1.cv);
    EXPECT_DOUBLE_EQ(state.cp, gold1.cp);
    EXPECT_NEAR(state.mu, gold1.mu, 2e-7);
    EXPECT_NEAR(state.k, gold1.k, 1e-3);
    EXPECT_DOUBLE_EQ(state.v, gold1.v);
    EXPECT_DOUBLE_EQ(state.s, gold1.s);
    EXPECT_DOUBLE_EQ(state.h, gold1.h);
    EXPECT_DOUBLE_EQ(state.w, gold1.w);
}

TEST(MethaneTest, rho_p)
{
    Methane fp;
    auto state = fp.rho_p(rho, p);

    EXPECT_DOUBLE_EQ(state.rho, gold1.rho);
    EXPECT_DOUBLE_EQ(state.T, gold1.T);
    EXPECT_DOUBLE_EQ(state.p, gold1.p);
    EXPECT_DOUBLE_EQ(state.u, gold1.u);
    EXPECT_DOUBLE_EQ(state.cv, gold1.cv);
    EXPECT_DOUBLE_EQ(state.cp, gold1.cp);
    EXPECT_NEAR(state.mu, gold1.mu, 2e-7);
    EXPECT_NEAR(state.k, gold1.k, 1e-3);
    EXPECT_DOUBLE_EQ(state.v, gold1.v);
    EXPECT_DOUBLE_EQ(state.s, gold1.s);
    EXPECT_DOUBLE_EQ(state.h, gold1.h);
    EXPECT_DOUBLE_EQ(state.w, gold1.w);
}

TEST(MethaneTest, p_T)
{
    Methane fp;
    auto state = fp.p_T(p, T);

    EXPECT_DOUBLE_EQ(state.rho, gold1.rho);
    EXPECT_DOUBLE_EQ(state.T, gold1.T);
    EXPECT_DOUBLE_EQ(state.p, gold1.p);
    EXPECT_DOUBLE_EQ(state.u, gold1.u);
    EXPECT_DOUBLE_EQ(state.cv, gold1.cv);
    EXPECT_DOUBLE_EQ(state.cp, gold1.cp);
    EXPECT_NEAR(state.mu, gold1.mu, 2e-7);
    EXPECT_NEAR(state.k, gold1.k, 1e-3);
    EXPECT_DOUBLE_EQ(state.v, gold1.v);
    EXPECT_DOUBLE_EQ(state.s, gold1.s);
    EXPECT_DOUBLE_EQ(state.h, gold1.h);
    EXPECT_DOUBLE_EQ(state.w, gold1.w);
}

TEST(MethaneTest, v_u)
{
    Methane fp;
    auto state0 = fp.p_T(p, T);
    State state = fp.v_u(state0.v, state0.u);

    EXPECT_DOUBLE_EQ(state.rho, gold1.rho);
    EXPECT_DOUBLE_EQ(state.T, gold1.T);
    EXPECT_DOUBLE_EQ(state.p, gold1.p);
    EXPECT_DOUBLE_EQ(state.u, gold1.u);
    EXPECT_DOUBLE_EQ(state.cv, gold1.cv);
    EXPECT_DOUBLE_EQ(state.cp, gold1.cp);
    EXPECT_NEAR(state.mu, gold1.mu, 2e-7);
    EXPECT_NEAR(state.k, gold1.k, 1e-3);
    EXPECT_DOUBLE_EQ(state.v, gold1.v);
    EXPECT_DOUBLE_EQ(state.s, gold1.s);
    EXPECT_DOUBLE_EQ(state.h, gold1.h);
    EXPECT_DOUBLE_EQ(state.w, gold1.w);
}
