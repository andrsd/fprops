#include "gtest/gtest.h"
#include <iomanip>
#include "fprops/nitrogen.h"

using namespace fprops;

namespace {

double rho = 12.074993451051515;
double v = 1. / rho;
double u = 205267.70993394594;
double T = 280.0;
double p = 1.0e6;
double cp = 1058.6154681901673;
double cv = 745.56949823705611;
double s = 6083.1854964583363;
double h = 288083.48983922758;
double c = 342.35848437431741;
double mu = 1.7090507109297636e-05;
double k = 0.024857419011067187;

// T = 280 K, p = 1 MPa
State gold1 = { u, v, rho, p, T, mu, cp, cv, s, k, h, c };

} // namespace

TEST(NitrogenTest, rho_T)
{
    Nitrogen fp;
    auto state = fp.rho_T(rho, T);

    EXPECT_DOUBLE_EQ(state.rho, gold1.rho);
    EXPECT_DOUBLE_EQ(state.T, gold1.T);
    EXPECT_DOUBLE_EQ(state.p, gold1.p);
    EXPECT_DOUBLE_EQ(state.u, gold1.u);
    EXPECT_DOUBLE_EQ(state.cv, gold1.cv);
    EXPECT_DOUBLE_EQ(state.cp, gold1.cp);
    EXPECT_DOUBLE_EQ(state.mu, gold1.mu);
    EXPECT_DOUBLE_EQ(state.k, gold1.k);
    EXPECT_NEAR(state.v, gold1.v, 1e-9);
    EXPECT_DOUBLE_EQ(state.s, gold1.s);
    EXPECT_DOUBLE_EQ(state.h, gold1.h);
    EXPECT_DOUBLE_EQ(state.w, gold1.w);
}

TEST(NitrogenTest, rho_p)
{
    Nitrogen fp;
    auto state = fp.rho_p(rho, p);

    EXPECT_DOUBLE_EQ(state.rho, gold1.rho);
    EXPECT_DOUBLE_EQ(state.T, gold1.T);
    EXPECT_DOUBLE_EQ(state.p, gold1.p);
    EXPECT_DOUBLE_EQ(state.u, gold1.u);
    EXPECT_DOUBLE_EQ(state.cv, gold1.cv);
    EXPECT_DOUBLE_EQ(state.cp, gold1.cp);
    EXPECT_DOUBLE_EQ(state.mu, gold1.mu);
    EXPECT_DOUBLE_EQ(state.k, gold1.k);
    EXPECT_NEAR(state.v, gold1.v, 1e-9);
    EXPECT_DOUBLE_EQ(state.s, gold1.s);
    EXPECT_DOUBLE_EQ(state.h, gold1.h);
    EXPECT_DOUBLE_EQ(state.w, gold1.w);
}

TEST(NitrogenTest, p_T)
{
    Nitrogen fp;
    auto state = fp.p_T(p, T);

    EXPECT_DOUBLE_EQ(state.rho, gold1.rho);
    EXPECT_DOUBLE_EQ(state.T, gold1.T);
    EXPECT_DOUBLE_EQ(state.p, gold1.p);
    EXPECT_DOUBLE_EQ(state.u, gold1.u);
    EXPECT_DOUBLE_EQ(state.cv, gold1.cv);
    EXPECT_DOUBLE_EQ(state.cp, gold1.cp);
    EXPECT_DOUBLE_EQ(state.mu, gold1.mu);
    EXPECT_DOUBLE_EQ(state.k, gold1.k);
    EXPECT_NEAR(state.v, gold1.v, 1e-9);
    EXPECT_DOUBLE_EQ(state.s, gold1.s);
    EXPECT_DOUBLE_EQ(state.h, gold1.h);
    EXPECT_DOUBLE_EQ(state.w, gold1.w);
}

TEST(NitrogenTest, v_u)
{
    Nitrogen fp;
    auto state = fp.v_u(v, u);

    EXPECT_NEAR(state.rho, gold1.rho, 1e-9);
    EXPECT_DOUBLE_EQ(state.T, gold1.T);
    EXPECT_NEAR(state.p, gold1.p, 1e-8);
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
