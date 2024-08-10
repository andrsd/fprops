#include "gtest/gtest.h"
#include "fprops/carbon_dioxide.h"

using namespace fprops;

namespace {

double p = 1e6;
double T = 280;
double rho = 20.19930913628047;
double v = 1. / rho;
double h = 480394.6994348873;
double u = 430888.0557544915;
double s = 2225.7418424216157;
double c = 252.3265654013213;
double cv = 670.919849374081;
double cp = 925.1800875340849;
double mu = 1.4126691599927136e-5;
double k = 0.01587080098501899;

// T = 280 K, p = 1 MPa
State gold1 = { u, v, rho, p, T, mu, cp, cv, s, k, h, c };

} // namespace

TEST(CarbonDioxide, rho_T)
{
    CarbonDioxide fp;
    auto state = fp.rho_T(rho, T);

    EXPECT_DOUBLE_EQ(state.rho, gold1.rho);
    EXPECT_DOUBLE_EQ(state.T, gold1.T);
    EXPECT_NEAR(state.p, gold1.p, 1e-2);
    EXPECT_DOUBLE_EQ(state.u, gold1.u);
    EXPECT_NEAR(state.cv, gold1.cv, 1e-5);
    EXPECT_NEAR(state.cp, gold1.cp, 2e-3);
    EXPECT_NEAR(state.mu, gold1.mu, 1e-7);
    EXPECT_NEAR(state.k, gold1.k, 3e-4);
    EXPECT_DOUBLE_EQ(state.v, gold1.v);
    EXPECT_DOUBLE_EQ(state.s, gold1.s);
    EXPECT_NEAR(state.h, gold1.h, 5e-3);
    EXPECT_NEAR(state.w, gold1.w, 1e-4);
}

TEST(CarbonDioxide, rho_p)
{
    CarbonDioxide fp;
    auto state = fp.rho_p(rho, p);

    EXPECT_DOUBLE_EQ(state.rho, gold1.rho);
    EXPECT_NEAR(state.T, gold1.T, 1e-5);
    EXPECT_NEAR(state.p, gold1.p, 1e-9);
    EXPECT_NEAR(state.u, gold1.u, 2e-3);
    EXPECT_NEAR(state.cv, gold1.cv, 1e-5);
    EXPECT_NEAR(state.cp, gold1.cp, 2e-3);
    EXPECT_NEAR(state.mu, gold1.mu, 1e-7);
    EXPECT_NEAR(state.k, gold1.k, 3e-4);
    EXPECT_DOUBLE_EQ(state.v, gold1.v);
    EXPECT_NEAR(state.s, gold1.s, 5e-6);
    EXPECT_NEAR(state.h, gold1.h, 5e-3);
    EXPECT_NEAR(state.w, gold1.w, 1e-4);
}

TEST(CarbonDioxide, p_T)
{
    CarbonDioxide fp;
    auto state = fp.p_T(p, T);

    EXPECT_NEAR(state.rho, gold1.rho, 2e-7);
    EXPECT_DOUBLE_EQ(state.T, gold1.T);
    EXPECT_NEAR(state.p, gold1.p, 1e-2);
    EXPECT_NEAR(state.u, gold1.u, 1e-4);
    EXPECT_NEAR(state.cv, gold1.cv, 1e-5);
    EXPECT_NEAR(state.cp, gold1.cp, 2e-3);
    EXPECT_NEAR(state.mu, gold1.mu, 1e-7);
    EXPECT_NEAR(state.k, gold1.k, 3e-4);
    EXPECT_NEAR(state.v, gold1.v, 1e-9);
    EXPECT_NEAR(state.s, gold1.s, 2e-6);
    EXPECT_NEAR(state.h, gold1.h, 5e-3);
    EXPECT_NEAR(state.w, gold1.w, 1e-4);
}

TEST(CarbonDioxide, v_u)
{
    CarbonDioxide fp;
    auto state = fp.v_u(v, u);

    EXPECT_DOUBLE_EQ(state.rho, gold1.rho);
    EXPECT_DOUBLE_EQ(state.T, gold1.T);
    EXPECT_NEAR(state.p, gold1.p, 1e-2);
    EXPECT_DOUBLE_EQ(state.u, gold1.u);
    EXPECT_NEAR(state.cv, gold1.cv, 1e-5);
    EXPECT_NEAR(state.cp, gold1.cp, 2e-3);
    EXPECT_NEAR(state.mu, gold1.mu, 1e-7);
    EXPECT_NEAR(state.k, gold1.k, 3e-4);
    EXPECT_DOUBLE_EQ(state.v, gold1.v);
    EXPECT_DOUBLE_EQ(state.s, gold1.s);
    EXPECT_NEAR(state.h, gold1.h, 5e-3);
    EXPECT_NEAR(state.w, gold1.w, 1e-4);
}
