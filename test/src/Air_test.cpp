#include "gtest/gtest.h"
#include "Air.h"

using namespace fprops;

TEST(Air, p_T)
{
    Air fp;

    double T = 300;
    double p = 101325;
    SinglePhaseFluidProperties::Props props = fp.p_T(p, T);

    EXPECT_DOUBLE_EQ(props.p, p);
    EXPECT_DOUBLE_EQ(props.T, T);
    EXPECT_DOUBLE_EQ(props.rho, 1.1769510785919943);
    EXPECT_DOUBLE_EQ(props.u, 3.400896182184841e5);
    EXPECT_DOUBLE_EQ(props.cv, 7.1756430588414378e2);
    EXPECT_DOUBLE_EQ(props.cp, 1.0058073121581035e3);
    EXPECT_DOUBLE_EQ(props.mu, 1.8537852519143559e-05);
    EXPECT_DOUBLE_EQ(props.k, 26.386197401795309e-3);
    EXPECT_DOUBLE_EQ(props.v, 0.84965298744304329);
    EXPECT_DOUBLE_EQ(props.s, 3.9227271289871092e3);
    EXPECT_DOUBLE_EQ(props.h, 4.2618070717115048e5);
    EXPECT_DOUBLE_EQ(props.w, 3.473335907030085e2);
}

TEST(Air, v_u)
{
    Air fp;

    double v = 0.84965298744304329;
    double u = 3.400896182184841e5;
    SinglePhaseFluidProperties::Props props = fp.v_u(v, u);

    EXPECT_DOUBLE_EQ(props.v, v);
    EXPECT_DOUBLE_EQ(props.u, u);
    EXPECT_DOUBLE_EQ(props.p, 101325);
    EXPECT_DOUBLE_EQ(props.T, 300);
    EXPECT_DOUBLE_EQ(props.rho, 1.1769510785919943);
    EXPECT_DOUBLE_EQ(props.cv, 7.1756430588414378e2);
    EXPECT_DOUBLE_EQ(props.cp, 1.0058073121581035e3);
    EXPECT_DOUBLE_EQ(props.mu, 1.8537852519143559e-05);
    EXPECT_DOUBLE_EQ(props.k, 26.386197401795309e-3);
    EXPECT_DOUBLE_EQ(props.s, 3.9227271289871092e3);
    EXPECT_DOUBLE_EQ(props.h, 4.2618070717115048e5);
    EXPECT_DOUBLE_EQ(props.w, 3.473335907030085e2);
}
