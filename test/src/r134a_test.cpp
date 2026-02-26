#include "gtest/gtest.h"
#include "fprops/r134a.h"

using namespace fprops;

// NOTE: The tolerance are quite loose. It could be because of:
// 1. The initial NL solve to get the (rho,T) point in Helmholz is not converged enough
// 2. The point for comparison is not picked in a reagion where we match better

namespace {

double rho = 4.083501;
double u = 409905.9099662705;
double v = 1. / rho;
double p = 101325;
double T = 310;
double mu = 1.2281512459243506e-05;
double cp = 869.2904624328216;
double cv = 779.7992677565012;
double s = 1934.1274797598596;
double k = 0.014337804704553059;
double h = 434719.1743801649;
double c = 164.81188238095876;

// T = 310 K, p = 101325 Pa
State gold1 = { u, v, rho, p, T, mu, cp, cv, s, k, h, c };

} // namespace

TEST(R134aTest, rho_T)
{
    R134a fp;
    auto state = fp.rho_T(rho, T);

    EXPECT_NEAR(state.rho, gold1.rho, 1e-10);
    EXPECT_NEAR(state.T, gold1.T, 1e-10);
    EXPECT_NEAR(state.p, gold1.p, 22);
    EXPECT_NEAR(state.u, gold1.u, 25);
    EXPECT_NEAR(state.cv, gold1.cv, 0.1);
    EXPECT_NEAR(state.cp, gold1.cp, 0.2);
    EXPECT_NEAR(state.mu, gold1.mu, 1e-7);
    EXPECT_NEAR(state.k, gold1.k, 5e-6);
    EXPECT_NEAR(state.v, gold1.v, 1e-10);
    EXPECT_NEAR(state.s, gold1.s, 1.);
    EXPECT_NEAR(state.h, gold1.h, 27);
    EXPECT_NEAR(state.w, gold1.w, 0.1);
}

TEST(R134aTest, rho_p)
{
    R134a fp;
    auto state = fp.rho_p(rho, p);

    EXPECT_NEAR(state.rho, gold1.rho, 1e-10);
    EXPECT_NEAR(state.T, gold1.T, 0.1);
    EXPECT_NEAR(state.p, gold1.p, 1e-10);
    EXPECT_NEAR(state.u, gold1.u, 10.);
    EXPECT_NEAR(state.cv, gold1.cv, 0.2);
    EXPECT_NEAR(state.cp, gold1.cp, 0.3);
    EXPECT_NEAR(state.mu, gold1.mu, 1e-7);
    EXPECT_NEAR(state.k, gold1.k, 5e-6);
    EXPECT_NEAR(state.v, gold1.v, 1e-10);
    EXPECT_NEAR(state.s, gold1.s, 0.5);
    EXPECT_NEAR(state.h, gold1.h, 10);
    EXPECT_NEAR(state.w, gold1.w, 0.01);
}

TEST(R134aTest, p_T)
{
    R134a fp;
    auto state = fp.p_T(p, T);

    EXPECT_NEAR(state.rho, gold1.rho, 1e-3);
    EXPECT_NEAR(state.T, gold1.T, 1e-10);
    EXPECT_NEAR(state.p, gold1.p, 1e-10);
    EXPECT_NEAR(state.u, gold1.u, 25);
    EXPECT_NEAR(state.cv, gold1.cv, 0.1);
    EXPECT_NEAR(state.cp, gold1.cp, 0.2);
    EXPECT_NEAR(state.mu, gold1.mu, 1e-7);
    EXPECT_NEAR(state.k, gold1.k, 5e-6);
    EXPECT_NEAR(state.v, gold1.v, 1e-4);
    EXPECT_NEAR(state.s, gold1.s, 0.6);
    EXPECT_NEAR(state.h, gold1.h, 30);
    EXPECT_NEAR(state.w, gold1.w, 0.05);
}

TEST(R134aTest, v_u)
{
    R134a fp;
    auto state = fp.v_u(v, u);

    EXPECT_NEAR(state.rho, gold1.rho, 1e-10);
    EXPECT_NEAR(state.T, gold1.T, 0.05);
    EXPECT_NEAR(state.p, gold1.p, 5.);
    EXPECT_NEAR(state.u, gold1.u, 1e-10);
    EXPECT_NEAR(state.cv, gold1.cv, 0.2);
    EXPECT_NEAR(state.cp, gold1.cp, 0.3);
    EXPECT_NEAR(state.mu, gold1.mu, 1e-7);
    EXPECT_NEAR(state.k, gold1.k, 5e-6);
    EXPECT_NEAR(state.v, gold1.v, 1e-10);
    EXPECT_NEAR(state.s, gold1.s, 1.);
    EXPECT_NEAR(state.h, gold1.h, 1.);
    EXPECT_NEAR(state.w, gold1.w, 0.01);
}
