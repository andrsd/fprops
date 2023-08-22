#include "gtest/gtest.h"
#include "Air.h"

using namespace fprops;

TEST(Air, rho_T)
{
    Air fp;

    double rho = 1.1769510785919943;
    double T = 300;
    SinglePhaseFluidProperties::Props props = fp.rho_T(rho, T);

    EXPECT_DOUBLE_EQ(props.rho, rho);
    EXPECT_DOUBLE_EQ(props.T, T);
    EXPECT_DOUBLE_EQ(props.p, 101325);
    EXPECT_DOUBLE_EQ(props.u, 340138.96058601438);
    EXPECT_DOUBLE_EQ(props.cv, 717.95390122400988);
    EXPECT_DOUBLE_EQ(props.cp, 1006.1967213037287);
    EXPECT_DOUBLE_EQ(props.mu, 1.8537852519143559e-05);
    EXPECT_DOUBLE_EQ(props.k, 26.386197401795309e-3);
    EXPECT_DOUBLE_EQ(props.v, 0.84965298744304329);
    EXPECT_DOUBLE_EQ(props.s, 3923.025553321696);
    EXPECT_DOUBLE_EQ(props.h, 426230.04953868076);
    EXPECT_DOUBLE_EQ(props.w, 347.30666259109006);
}

TEST(Air, p_T)
{
    Air fp;

    double T = 300;
    double p = 101325;
    SinglePhaseFluidProperties::Props props = fp.p_T(p, T);

    EXPECT_DOUBLE_EQ(props.p, p);
    EXPECT_DOUBLE_EQ(props.T, T);
    EXPECT_DOUBLE_EQ(props.rho, 1.1769510785919943);
    EXPECT_DOUBLE_EQ(props.u, 340138.96058601438);
    EXPECT_DOUBLE_EQ(props.cv, 717.95390122400988);
    EXPECT_DOUBLE_EQ(props.cp, 1006.1967213037287);
    EXPECT_DOUBLE_EQ(props.mu, 1.8537852519143559e-05);
    EXPECT_DOUBLE_EQ(props.k, 26.386197401795309e-3);
    EXPECT_DOUBLE_EQ(props.v, 0.84965298744304329);
    EXPECT_DOUBLE_EQ(props.s, 3923.025553321696);
    EXPECT_DOUBLE_EQ(props.h, 426230.04953868076);
    EXPECT_DOUBLE_EQ(props.w, 347.30666259109006);
}

TEST(Air, v_u)
{
    Air fp;

    double v = 0.84965298744304329;
    double u = 340138.96058601438;
    SinglePhaseFluidProperties::Props props = fp.v_u(v, u);

    EXPECT_DOUBLE_EQ(props.v, v);
    EXPECT_DOUBLE_EQ(props.u, u);
    EXPECT_DOUBLE_EQ(props.p, 101325);
    EXPECT_DOUBLE_EQ(props.T, 300);
    EXPECT_DOUBLE_EQ(props.rho, 1.1769510785919943);
    EXPECT_DOUBLE_EQ(props.cv, 717.95390122400988);
    EXPECT_DOUBLE_EQ(props.cp, 1006.1967213037287);
    EXPECT_DOUBLE_EQ(props.mu, 1.8537852519143559e-05);
    EXPECT_DOUBLE_EQ(props.k, 26.386197401795309e-3);
    EXPECT_DOUBLE_EQ(props.s, 3923.025553321696);
    EXPECT_DOUBLE_EQ(props.h, 426230.04953868076);
    EXPECT_DOUBLE_EQ(props.w, 347.30666259109006);
}
