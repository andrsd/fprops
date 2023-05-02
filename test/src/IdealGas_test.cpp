#include "gtest/gtest.h"
#include "IdealGas.h"

using namespace fprops;

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

    EXPECT_DOUBLE_EQ(props.p, p);
    EXPECT_DOUBLE_EQ(props.T, T);
    EXPECT_DOUBLE_EQ(props.rho, 0.89892258591830565);
    EXPECT_DOUBLE_EQ(props.u, 2.8179567848017247e5);
    EXPECT_DOUBLE_EQ(props.cv, 7.16763775862069e2);
    EXPECT_DOUBLE_EQ(props.cp, 1.0034692862068968e3);
    EXPECT_DOUBLE_EQ(props.mu, 18.23e-6);
    EXPECT_DOUBLE_EQ(props.k, 25.68e-3);
    EXPECT_DOUBLE_EQ(props.v, 1.1124428462084279);
    EXPECT_DOUBLE_EQ(props.s, 2.6903243258630837e3);
    EXPECT_DOUBLE_EQ(props.h, 3.9451394987224141e5);
    EXPECT_DOUBLE_EQ(props.w, 3.9724750464779078e2);
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

    EXPECT_DOUBLE_EQ(props.v, v);
    EXPECT_DOUBLE_EQ(props.u, u);
    EXPECT_DOUBLE_EQ(props.p, 101325);
    EXPECT_DOUBLE_EQ(props.T, 393.15);
    EXPECT_DOUBLE_EQ(props.rho, 0.89892258591830565);
    EXPECT_DOUBLE_EQ(props.cv, 7.16763775862069e2);
    EXPECT_DOUBLE_EQ(props.cp, 1.0034692862068968e3);
    EXPECT_DOUBLE_EQ(props.mu, 18.23e-6);
    EXPECT_DOUBLE_EQ(props.k, 25.68e-3);
    EXPECT_DOUBLE_EQ(props.s, 2.6903243258630837e3);
    EXPECT_DOUBLE_EQ(props.h, 3.9451394987224141e5);
    EXPECT_DOUBLE_EQ(props.w, 3.9724750464779078e2);
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
