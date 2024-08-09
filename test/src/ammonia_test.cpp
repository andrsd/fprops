#include "gtest/gtest.h"
#include "fprops/ammonia.h"

using namespace fprops;

namespace {

double u = 1505100.2440409528;
double rho = 7.699184746208876;
double v = 1. / rho;
double p = 1.e6;
double T = 300.;
double mu = 9.911559307081025e-06;
double cp = 3097.3320775845727;
double cv = 2121.9092763848703;
double s = 5823.126832132577;
double k = 0.026307739861873722;
double h = 1634984.1256418484;
double w = 406.68113581057514;

State gold1 = { u, v, rho, p, T, mu, cp, cv, s, k, h, w };

} // namespace

TEST(Ammonia, rho_T)
{
    Ammonia fp;
    auto state = fp.rho_T(rho, T);

    EXPECT_NEAR(state.rho, gold1.rho, 1e-7);
    EXPECT_NEAR(state.T, gold1.T, 1e-7);
    EXPECT_NEAR(state.p, gold1.p, 1e-7);
    EXPECT_NEAR(state.u, gold1.u, 1e-3);
    EXPECT_NEAR(state.cv, gold1.cv, 3e-6);
    EXPECT_NEAR(state.cp, gold1.cp, 3e-6);
    EXPECT_NEAR(state.mu, gold1.mu, 5e-9);
    EXPECT_NEAR(state.k, gold1.k, 1e-7);
    EXPECT_NEAR(state.v, gold1.v, 1e-7);
    EXPECT_NEAR(state.s, gold1.s, 1e-6);
    EXPECT_NEAR(state.h, gold1.h, 3e-4);
    EXPECT_NEAR(state.w, gold1.w, 1e-7);
}

TEST(Ammonia, rho_p)
{
    Ammonia fp;
    auto state = fp.rho_p(rho, p);

    EXPECT_NEAR(state.rho, gold1.rho, 1e-7);
    EXPECT_NEAR(state.T, gold1.T, 1e-7);
    EXPECT_NEAR(state.p, gold1.p, 1e-7);
    EXPECT_NEAR(state.u, gold1.u, 1e-3);
    EXPECT_NEAR(state.cv, gold1.cv, 3e-6);
    EXPECT_NEAR(state.cp, gold1.cp, 3e-6);
    EXPECT_NEAR(state.mu, gold1.mu, 5e-9);
    EXPECT_NEAR(state.k, gold1.k, 1e-7);
    EXPECT_NEAR(state.v, gold1.v, 1e-7);
    EXPECT_NEAR(state.s, gold1.s, 1e-6);
    EXPECT_NEAR(state.h, gold1.h, 3e-4);
    EXPECT_NEAR(state.w, gold1.w, 1e-7);
}

TEST(Ammonia, p_T)
{
    Ammonia fp;
    auto state = fp.p_T(p, T);

    EXPECT_NEAR(state.rho, gold1.rho, 1e-6);
    EXPECT_NEAR(state.T, gold1.T, 1e-6);
    EXPECT_NEAR(state.p, gold1.p, 1e-6);
    EXPECT_NEAR(state.u, gold1.u, 1e-3);
    EXPECT_NEAR(state.cv, gold1.cv, 3e-6);
    EXPECT_NEAR(state.cp, gold1.cp, 3e-6);
    EXPECT_NEAR(state.mu, gold1.mu, 5e-9);
    EXPECT_NEAR(state.k, gold1.k, 1e-6);
    EXPECT_NEAR(state.v, gold1.v, 1e-7);
    EXPECT_NEAR(state.s, gold1.s, 1e-6);
    EXPECT_NEAR(state.h, gold1.h, 3e-4);
    EXPECT_NEAR(state.w, gold1.w, 1e-7);
}

TEST(Ammonia, v_u)
{
    Ammonia fp;
    auto state = fp.v_u(v, u);

    EXPECT_NEAR(state.rho, gold1.rho, 1e-7);
    EXPECT_NEAR(state.T, gold1.T, 1e-7);
    EXPECT_NEAR(state.p, gold1.p, 1e-3);
    EXPECT_NEAR(state.u, gold1.u, 1e-7);
    EXPECT_NEAR(state.cv, gold1.cv, 2e-6);
    EXPECT_NEAR(state.cp, gold1.cp, 1e-6);
    EXPECT_NEAR(state.mu, gold1.mu, 1e-7);
    EXPECT_NEAR(state.k, gold1.k, 1e-7);
    EXPECT_NEAR(state.v, gold1.v, 1e-7);
    EXPECT_NEAR(state.s, gold1.s, 3e-7);
    EXPECT_NEAR(state.h, gold1.h, 3e-4);
    EXPECT_NEAR(state.w, gold1.w, 1e-7);
}
