#include "gtest/gtest.h"
#include "IdealGas.h"

using namespace fprops;

namespace {

SinglePhaseFluidProperties::Props gold1 = { 2.8179567848017247e5,
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

TEST(IdealGas, rho_T)
{
    double gamma = 1.4;
    double molar_mass = 29.0e-3;
    IdealGas fp(gamma, molar_mass);
    fp.set_mu(18.23e-6);
    fp.set_k(25.68e-3);

    double rho = 0.89892258591830565;
    double T = 120. + 273.15;
    SinglePhaseFluidProperties::Props props = fp.rho_T(rho, T);

    EXPECT_DOUBLE_EQ(props.rho, gold1.rho);
    EXPECT_DOUBLE_EQ(props.T, gold1.T);
    EXPECT_DOUBLE_EQ(props.p, gold1.p);
    EXPECT_DOUBLE_EQ(props.u, gold1.u);
    EXPECT_DOUBLE_EQ(props.cv, gold1.cv);
    EXPECT_DOUBLE_EQ(props.cp, gold1.cp);
    EXPECT_DOUBLE_EQ(props.mu, gold1.mu);
    EXPECT_DOUBLE_EQ(props.k, gold1.k);
    EXPECT_DOUBLE_EQ(props.v, gold1.v);
    EXPECT_DOUBLE_EQ(props.s, gold1.s);
    EXPECT_DOUBLE_EQ(props.h, gold1.h);
    EXPECT_DOUBLE_EQ(props.w, gold1.w);
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
    SinglePhaseFluidProperties::Props props = fp.rho_p(rho, p);

    EXPECT_DOUBLE_EQ(props.rho, gold1.rho);
    EXPECT_DOUBLE_EQ(props.T, gold1.T);
    EXPECT_DOUBLE_EQ(props.p, gold1.p);
    EXPECT_DOUBLE_EQ(props.u, gold1.u);
    EXPECT_DOUBLE_EQ(props.cv, gold1.cv);
    EXPECT_DOUBLE_EQ(props.cp, gold1.cp);
    EXPECT_DOUBLE_EQ(props.mu, gold1.mu);
    EXPECT_DOUBLE_EQ(props.k, gold1.k);
    EXPECT_DOUBLE_EQ(props.v, gold1.v);
    EXPECT_DOUBLE_EQ(props.s, gold1.s);
    EXPECT_DOUBLE_EQ(props.h, gold1.h);
    EXPECT_DOUBLE_EQ(props.w, gold1.w);
}

TEST(IdealGas, rho_p_incorrect)
{
    double gamma = 1.4;
    double molar_mass = 29.0e-3;
    IdealGas fp(gamma, molar_mass);

    EXPECT_THROW(auto p = fp.rho_p(-1, 1e5), std::domain_error);
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
    SinglePhaseFluidProperties::Props props = fp.p_T(p, T);

    EXPECT_DOUBLE_EQ(props.rho, gold1.rho);
    EXPECT_DOUBLE_EQ(props.T, gold1.T);
    EXPECT_DOUBLE_EQ(props.p, gold1.p);
    EXPECT_DOUBLE_EQ(props.u, gold1.u);
    EXPECT_DOUBLE_EQ(props.cv, gold1.cv);
    EXPECT_DOUBLE_EQ(props.cp, gold1.cp);
    EXPECT_DOUBLE_EQ(props.mu, gold1.mu);
    EXPECT_DOUBLE_EQ(props.k, gold1.k);
    EXPECT_DOUBLE_EQ(props.v, gold1.v);
    EXPECT_DOUBLE_EQ(props.s, gold1.s);
    EXPECT_DOUBLE_EQ(props.h, gold1.h);
    EXPECT_DOUBLE_EQ(props.w, gold1.w);
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
    SinglePhaseFluidProperties::Props props = fp.v_u(v, u);

    EXPECT_DOUBLE_EQ(props.rho, gold1.rho);
    EXPECT_DOUBLE_EQ(props.T, gold1.T);
    EXPECT_DOUBLE_EQ(props.p, gold1.p);
    EXPECT_DOUBLE_EQ(props.u, gold1.u);
    EXPECT_DOUBLE_EQ(props.cv, gold1.cv);
    EXPECT_DOUBLE_EQ(props.cp, gold1.cp);
    EXPECT_DOUBLE_EQ(props.mu, gold1.mu);
    EXPECT_DOUBLE_EQ(props.k, gold1.k);
    EXPECT_DOUBLE_EQ(props.v, gold1.v);
    EXPECT_DOUBLE_EQ(props.s, gold1.s);
    EXPECT_DOUBLE_EQ(props.h, gold1.h);
    EXPECT_DOUBLE_EQ(props.w, gold1.w);
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
    SinglePhaseFluidProperties::Props props = fp.h_s(h, s);

    EXPECT_DOUBLE_EQ(props.rho, gold1.rho);
    EXPECT_NEAR(props.T, gold1.T, 1e-12);
    EXPECT_NEAR(props.p, gold1.p, 1e-9);
    EXPECT_NEAR(props.u, gold1.u, 1e-9);
    EXPECT_DOUBLE_EQ(props.cv, gold1.cv);
    EXPECT_DOUBLE_EQ(props.cp, gold1.cp);
    EXPECT_DOUBLE_EQ(props.mu, gold1.mu);
    EXPECT_DOUBLE_EQ(props.k, gold1.k);
    EXPECT_DOUBLE_EQ(props.v, gold1.v);
    EXPECT_DOUBLE_EQ(props.s, gold1.s);
    EXPECT_DOUBLE_EQ(props.h, gold1.h);
    EXPECT_NEAR(props.w, gold1.w, 1e-9);
}

TEST(IdealGas, v_u_incorrect)
{
    double gamma = 1.4;
    double molar_mass = 29.0e-3;
    IdealGas fp(gamma, molar_mass);

    EXPECT_THROW(auto p = fp.v_u(-1, 1), std::domain_error);
    EXPECT_THROW(auto p = fp.v_u(1, -1), std::domain_error);
}

TEST(IdealGas, p_T_incorrect)
{
    double gamma = 1.4;
    double molar_mass = 29.0e-3;
    IdealGas fp(gamma, molar_mass);

    EXPECT_THROW(auto p = fp.p_T(1e5, -1), std::domain_error);
}
