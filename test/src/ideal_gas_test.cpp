#include "gtest/gtest.h"
#include "exception_test_macros.h"
#include "fprops/ideal_gas.h"

using namespace fprops;

namespace {

State gold1 = { 2.8179567848017247e5,
                1.1124428462084279,
                0.89892258591830565,
                101325,
                393.15,
                18.23e-6,
                1.0034692862068968e3,
                7.16763775862069e2,
                2.6903243258630837e3,
                25.68e-3,
                3.9451394987224141e5,
                3.9724750464779078e2 };

}

TEST(IdealGas, api)
{
    double gamma = 1.4;
    double molar_mass = 29.0e-3;
    double mu = 18.23e-6;
    double k = 25.68e-3;
    IdealGas fp(gamma, molar_mass);
    fp.set_mu(mu);
    fp.set_k(k);

    EXPECT_DOUBLE_EQ(fp.gamma(), gamma);
    EXPECT_NEAR(fp.R_specific(), 286.7055103448, 1e-10);
    EXPECT_DOUBLE_EQ(fp.molar_mass(), molar_mass);
    EXPECT_NEAR(fp.cp(), 1003.4692862068, 1e-10);
    EXPECT_NEAR(fp.cv(), 716.7637758620, 1e-10);
    EXPECT_DOUBLE_EQ(fp.mu(), mu);
    EXPECT_DOUBLE_EQ(fp.k(), k);
}

TEST(IdealGas, rho_T)
{
    double gamma = 1.4;
    double molar_mass = 29.0e-3;
    IdealGas fp(gamma, molar_mass);
    fp.set_mu(18.23e-6);
    fp.set_k(25.68e-3);

    double rho = 0.89892258591830565;
    double T = 120. + 273.15;
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

TEST(IdealGas, rho_p)
{
    double gamma = 1.4;
    double molar_mass = 29.0e-3;
    IdealGas fp(gamma, molar_mass);
    fp.set_mu(18.23e-6);
    fp.set_k(25.68e-3);

    double rho = 0.89892258591830565;
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

TEST(IdealGas, rho_p_incorrect)
{
    double gamma = 1.4;
    double molar_mass = 29.0e-3;
    IdealGas fp(gamma, molar_mass);

    EXPECT_THROW_MSG(auto st = fp.rho_p(-1, 1e5), "Negative density");
}

TEST(IdealGas, p_T)
{
    double gamma = 1.4;
    double molar_mass = 29.0e-3;
    IdealGas fp(gamma, molar_mass);
    fp.set_mu(18.23e-6);
    fp.set_k(25.68e-3);

    double T = 120. + 273.15;
    double p = 101325;
    State state = fp.p_T(p, T);

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

TEST(IdealGas, v_u)
{
    double gamma = 1.4;
    double molar_mass = 29.0e-3;
    IdealGas fp(gamma, molar_mass);
    fp.set_mu(18.23e-6);
    fp.set_k(25.68e-3);

    double v = 1.1124428462084279;
    double u = 2.8179567848017247e5;
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

TEST(IdealGas, h_s)
{
    double gamma = 1.4;
    double molar_mass = 29.0e-3;
    IdealGas fp(gamma, molar_mass);
    fp.set_mu(18.23e-6);
    fp.set_k(25.68e-3);

    double h = 3.9451394987224141e5;
    double s = 2.6903243258630837e3;
    auto state = fp.h_s(h, s);

    EXPECT_DOUBLE_EQ(state.rho, gold1.rho);
    EXPECT_NEAR(state.T, gold1.T, 1e-12);
    EXPECT_NEAR(state.p, gold1.p, 1e-9);
    EXPECT_NEAR(state.u, gold1.u, 1e-9);
    EXPECT_DOUBLE_EQ(state.cv, gold1.cv);
    EXPECT_DOUBLE_EQ(state.cp, gold1.cp);
    EXPECT_DOUBLE_EQ(state.mu, gold1.mu);
    EXPECT_DOUBLE_EQ(state.k, gold1.k);
    EXPECT_DOUBLE_EQ(state.v, gold1.v);
    EXPECT_DOUBLE_EQ(state.s, gold1.s);
    EXPECT_DOUBLE_EQ(state.h, gold1.h);
    EXPECT_NEAR(state.w, gold1.w, 1e-9);
}

TEST(IdealGas, v_h)
{
    double gamma = 1.4;
    double molar_mass = 29.0e-3;
    IdealGas fp(gamma, molar_mass);
    fp.set_mu(18.23e-6);
    fp.set_k(25.68e-3);

    double v = 1.1124428462084279;
    double h = 3.9451394987224141e5;
    auto state = fp.v_h(v, h);

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

TEST(IdealGas, v_u_incorrect)
{
    double gamma = 1.4;
    double molar_mass = 29.0e-3;
    IdealGas fp(gamma, molar_mass);

    EXPECT_THROW_MSG(auto st = fp.v_u(-1, 1), "Negative specific volume");
    EXPECT_THROW_MSG(auto st = fp.v_u(1, -1), "Negative internal energy");
}

TEST(IdealGas, p_T_incorrect)
{
    double gamma = 1.4;
    double molar_mass = 29.0e-3;
    IdealGas fp(gamma, molar_mass);

    EXPECT_THROW_MSG(auto st = fp.p_T(1e5, -1), "Negative temperature");
}
TEST(IdealGas, v_h_incorrect)
{
    double gamma = 1.4;
    double molar_mass = 29.0e-3;
    IdealGas fp(gamma, molar_mass);

    EXPECT_THROW_MSG(auto st = fp.v_h(-1, 1), "Negative specific volume");
}
