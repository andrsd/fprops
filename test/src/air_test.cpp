#include "gtest/gtest.h"
#include "fprops/air.h"

using namespace fprops;

namespace {

// T = 300 K, p = 101325 Pa
State gold1 = { 340138.96058601438,
                0.84965298744304329,
                1.1769510785919943,
                101325,
                300,
                1.8537852519143559e-05,
                1006.1967213037287,
                717.95390122400988,
                3923.025553321696,
                26.386197401795309e-3,
                426230.04953868076,
                347.30666259109006 };

} // namespace

TEST(Air, rho_T)
{
    Air fp;

    double rho = 1.1769510785919943;
    double T = 300;
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

TEST(Air, rho_p)
{
    Air fp;

    double rho = 1.1769510785919943;
    double p = 101325;
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

TEST(Air, p_T)
{
    Air fp;

    double T = 300;
    double p = 101325;
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

TEST(Air, v_u)
{
    Air fp;

    double v = 0.84965298744304329;
    double u = 340138.96058601438;
    auto state = fp.v_u(v, u);

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
