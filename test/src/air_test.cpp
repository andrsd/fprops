#include "gtest/gtest.h"
#include "fprops/air.h"

using namespace fprops;

namespace {

double rho = 1.1769510785919943;
double u = 340138.96058601438;
double v = 1. / rho;
double p = 101325;
double T = 300;
double mu = 1.8537852519143559e-05;
double cp = 1006.1967213037287;
double cv = 717.95390122400988;
double s = 3923.025553321696;
double k = 26.386197401795309e-3;
double h = 426230.04953868076;
double c = 347.30666259109006;

// T = 300 K, p = 101325 Pa
State gold1 = { u, v, rho, p, T, mu, cp, cv, s, k, h, c };

} // namespace

TEST(Air, rho_T)
{
    Air fp;
    auto state = fp.rho_T(rho, T);

    EXPECT_DOUBLE_EQ(state.rho, gold1.rho);
    EXPECT_DOUBLE_EQ(state.T, gold1.T);
    EXPECT_DOUBLE_EQ(state.p, gold1.p);
    EXPECT_DOUBLE_EQ(state.u, gold1.u);
    EXPECT_DOUBLE_EQ(state.cv, gold1.cv);
    EXPECT_DOUBLE_EQ(state.cp, gold1.cp);
    EXPECT_NEAR(state.mu, gold1.mu, 1e-8);
    EXPECT_NEAR(state.k, gold1.k, 5e-6);
    EXPECT_DOUBLE_EQ(state.v, gold1.v);
    EXPECT_DOUBLE_EQ(state.s, gold1.s);
    EXPECT_DOUBLE_EQ(state.h, gold1.h);
    EXPECT_DOUBLE_EQ(state.w, gold1.w);
}

TEST(Air, rho_p)
{
    Air fp;
    auto state = fp.rho_p(rho, p);

    EXPECT_DOUBLE_EQ(state.rho, gold1.rho);
    EXPECT_DOUBLE_EQ(state.T, gold1.T);
    EXPECT_DOUBLE_EQ(state.p, gold1.p);
    EXPECT_DOUBLE_EQ(state.u, gold1.u);
    EXPECT_DOUBLE_EQ(state.cv, gold1.cv);
    EXPECT_DOUBLE_EQ(state.cp, gold1.cp);
    EXPECT_NEAR(state.mu, gold1.mu, 1e-8);
    EXPECT_NEAR(state.k, gold1.k, 5e-6);
    EXPECT_DOUBLE_EQ(state.v, gold1.v);
    EXPECT_DOUBLE_EQ(state.s, gold1.s);
    EXPECT_DOUBLE_EQ(state.h, gold1.h);
    EXPECT_DOUBLE_EQ(state.w, gold1.w);
}

TEST(Air, p_T)
{
    Air fp;
    auto state = fp.p_T(p, T);

    EXPECT_DOUBLE_EQ(state.rho, gold1.rho);
    EXPECT_DOUBLE_EQ(state.T, gold1.T);
    EXPECT_DOUBLE_EQ(state.p, gold1.p);
    EXPECT_DOUBLE_EQ(state.u, gold1.u);
    EXPECT_DOUBLE_EQ(state.cv, gold1.cv);
    EXPECT_DOUBLE_EQ(state.cp, gold1.cp);
    EXPECT_NEAR(state.mu, gold1.mu, 1e-8);
    EXPECT_NEAR(state.k, gold1.k, 5e-6);
    EXPECT_DOUBLE_EQ(state.v, gold1.v);
    EXPECT_DOUBLE_EQ(state.s, gold1.s);
    EXPECT_DOUBLE_EQ(state.h, gold1.h);
    EXPECT_DOUBLE_EQ(state.w, gold1.w);
}

TEST(Air, v_u)
{
    Air fp;
    auto state = fp.v_u(v, u);

    EXPECT_DOUBLE_EQ(state.rho, gold1.rho);
    EXPECT_DOUBLE_EQ(state.T, gold1.T);
    EXPECT_DOUBLE_EQ(state.p, gold1.p);
    EXPECT_DOUBLE_EQ(state.u, gold1.u);
    EXPECT_DOUBLE_EQ(state.cv, gold1.cv);
    EXPECT_DOUBLE_EQ(state.cp, gold1.cp);
    EXPECT_NEAR(state.mu, gold1.mu, 1e-8);
    EXPECT_NEAR(state.k, gold1.k, 5e-6);
    EXPECT_DOUBLE_EQ(state.v, gold1.v);
    EXPECT_DOUBLE_EQ(state.s, gold1.s);
    EXPECT_DOUBLE_EQ(state.h, gold1.h);
    EXPECT_DOUBLE_EQ(state.w, gold1.w);
}
