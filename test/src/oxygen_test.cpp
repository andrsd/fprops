#include "gtest/gtest.h"
#include "fprops/oxygen.h"

using namespace fprops;

namespace {

double rho = 1.3006898082601743;
double v = 1. / rho;
double u = 194810.52158664705;
double T = 300.;
double p = 101325.;
double cp = 919.887200720042;
double cv = 658.7297701599554;
double s = 6412.852423886543;
double h = 272711.49332989;
double c = 329.72289112160996;
double mu = 2.0652424277518874e-5;
double k = 0.02648596367516763;

// T = 300 K, p = 101325 Pa
State gold = { u, v, rho, p, T, mu, cp, cv, s, k, h, c };

} // namespace

TEST(Oxygen, rho_T)
{
    Oxygen fp;
    auto state = fp.rho_T(rho, T);

    EXPECT_DOUBLE_EQ(state.rho, gold.rho);
    EXPECT_DOUBLE_EQ(state.T, gold.T);
    EXPECT_DOUBLE_EQ(state.p, gold.p);
    EXPECT_DOUBLE_EQ(state.u, gold.u);
    EXPECT_DOUBLE_EQ(state.cv, gold.cv);
    EXPECT_DOUBLE_EQ(state.cp, gold.cp);
    EXPECT_DOUBLE_EQ(state.mu, gold.mu);
    EXPECT_NEAR(state.k, gold.k, 1e-8);
    EXPECT_DOUBLE_EQ(state.v, gold.v);
    EXPECT_DOUBLE_EQ(state.s, gold.s);
    EXPECT_DOUBLE_EQ(state.h, gold.h);
    EXPECT_DOUBLE_EQ(state.w, gold.w);
}

TEST(Oxygen, rho_p)
{
    Oxygen fp;
    auto state = fp.rho_p(rho, p);

    EXPECT_DOUBLE_EQ(state.rho, gold.rho);
    EXPECT_DOUBLE_EQ(state.T, gold.T);
    EXPECT_DOUBLE_EQ(state.p, gold.p);
    EXPECT_DOUBLE_EQ(state.u, gold.u);
    EXPECT_DOUBLE_EQ(state.cv, gold.cv);
    EXPECT_DOUBLE_EQ(state.cp, gold.cp);
    EXPECT_DOUBLE_EQ(state.mu, gold.mu);
    EXPECT_NEAR(state.k, gold.k, 1e-8);
    EXPECT_DOUBLE_EQ(state.v, gold.v);
    EXPECT_DOUBLE_EQ(state.s, gold.s);
    EXPECT_DOUBLE_EQ(state.h, gold.h);
    EXPECT_DOUBLE_EQ(state.w, gold.w);
}

TEST(Oxygen, p_T)
{
    Oxygen fp;
    auto state = fp.p_T(p, T);

    EXPECT_DOUBLE_EQ(state.rho, gold.rho);
    EXPECT_DOUBLE_EQ(state.T, gold.T);
    EXPECT_DOUBLE_EQ(state.p, gold.p);
    EXPECT_DOUBLE_EQ(state.u, gold.u);
    EXPECT_DOUBLE_EQ(state.cv, gold.cv);
    EXPECT_DOUBLE_EQ(state.cp, gold.cp);
    EXPECT_DOUBLE_EQ(state.mu, gold.mu);
    EXPECT_NEAR(state.k, gold.k, 1e-8);
    EXPECT_DOUBLE_EQ(state.v, gold.v);
    EXPECT_DOUBLE_EQ(state.s, gold.s);
    EXPECT_DOUBLE_EQ(state.h, gold.h);
    EXPECT_DOUBLE_EQ(state.w, gold.w);
}

TEST(Oxygen, v_u)
{
    Oxygen fp;
    auto state = fp.v_u(v, u);

    EXPECT_DOUBLE_EQ(state.rho, gold.rho);
    EXPECT_DOUBLE_EQ(state.T, gold.T);
    EXPECT_DOUBLE_EQ(state.p, gold.p);
    EXPECT_DOUBLE_EQ(state.u, gold.u);
    EXPECT_DOUBLE_EQ(state.cv, gold.cv);
    EXPECT_DOUBLE_EQ(state.cp, gold.cp);
    EXPECT_DOUBLE_EQ(state.mu, gold.mu);
    EXPECT_NEAR(state.k, gold.k, 1e-8);
    EXPECT_DOUBLE_EQ(state.v, gold.v);
    EXPECT_DOUBLE_EQ(state.s, gold.s);
    EXPECT_DOUBLE_EQ(state.h, gold.h);
    EXPECT_DOUBLE_EQ(state.w, gold.w);
}
