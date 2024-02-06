#include "gmock/gmock.h"
#include "fprops/InterpolatedFluidProperties.h"
#include <filesystem>

using namespace fprops;
using namespace std::filesystem;

// T = 280 K, p = 1 MPa
State gold1 = { 205267.70993394594,
                0.082815779905282,
                12.074993451051515,
                1.0e6,
                280.0,
                1.7090507109297636e-05,
                1058.6154681901673,
                745.56949823705611,
                6083.1854964583363,
                0.024857419011067187,
                288083.48983922758,
                342.35848437431741 };

TEST(InterpolatedFluidPropertiesTest, p_T)
{
    auto file_name = path(FPROPS_UNIT_TESTS_ROOT) / path("assets") / path("n2-smol.fprops.h5");
    InterpolatedFluidProperties fp;
    fp.load(file_name.string());

    double T = 280.0;
    double p = 1.0e6;
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

TEST(InterpolatedFluidPropertiesTest, rho_T)
{
    auto file_name = path(FPROPS_UNIT_TESTS_ROOT) / path("assets") / path("n2-smol.fprops.h5");
    InterpolatedFluidProperties fp;
    fp.load(file_name.string());

    double rho = 12.074993451051515;
    double T = 280.0;
    auto state = fp.rho_T(rho, T);

    EXPECT_NEAR(state.rho, gold1.rho, 1e-5);
    EXPECT_NEAR(state.T, gold1.T, 1e-9);
    EXPECT_NEAR(state.p, gold1.p, 4);
    EXPECT_NEAR(state.u, gold1.u, 1.25);
    EXPECT_NEAR(state.cv, gold1.cv, 1e-5);
    EXPECT_NEAR(state.cp, gold1.cp, 1e-4);
    EXPECT_NEAR(state.mu, gold1.mu, 1e-6);
    EXPECT_NEAR(state.k, gold1.k, 1e-6);
    EXPECT_NEAR(state.v, gold1.v, 1e-7);
    EXPECT_NEAR(state.s, gold1.s, 0.005);
    EXPECT_NEAR(state.h, gold1.h, 1.8);
    EXPECT_NEAR(state.w, gold1.w, 0.002);
}

TEST(InterpolatedFluidPropertiesTest, rho_p)
{
    auto file_name = path(FPROPS_UNIT_TESTS_ROOT) / path("assets") / path("n2-smol.fprops.h5");
    InterpolatedFluidProperties fp;
    fp.load(file_name.string());

    double rho = 12.074993451051515;
    double p = 1.0e6;
    auto state = fp.rho_p(rho, p);

    EXPECT_NEAR(state.rho, gold1.rho, 0.);
    EXPECT_NEAR(state.T, gold1.T, 0.012);
    EXPECT_NEAR(state.p, gold1.p, 0.);
    EXPECT_NEAR(state.u, gold1.u, 9.);
    EXPECT_NEAR(state.cv, gold1.cv, 0.0003);
    EXPECT_NEAR(state.cp, gold1.cp, 1e-4);
    EXPECT_NEAR(state.mu, gold1.mu, 1e-9);
    EXPECT_NEAR(state.k, gold1.k, 1e-6);
    EXPECT_NEAR(state.v, gold1.v, 1e-7);
    EXPECT_NEAR(state.s, gold1.s, 0.03);
    EXPECT_NEAR(state.h, gold1.h, 12.6);
    EXPECT_NEAR(state.w, gold1.w, 0.007);
}

TEST(InterpolatedFluidPropertiesTest, v_u)
{
    auto file_name = path(FPROPS_UNIT_TESTS_ROOT) / path("assets") / path("n2-smol.fprops.h5");
    InterpolatedFluidProperties fp;
    fp.load(file_name.string());

    double v = 0.082815779905282;
    double u = 205267.70993394594;
    auto state = fp.v_u(v, u);

    EXPECT_NEAR(state.rho, gold1.rho, 0.16);
    EXPECT_NEAR(state.T, gold1.T, 0.16);
    EXPECT_NEAR(state.p, gold1.p, 13500);
    EXPECT_NEAR(state.u, gold1.u, 1e-9);
    EXPECT_NEAR(state.cv, gold1.cv, 0.04);
    EXPECT_NEAR(state.cp, gold1.cp, 0.24);
    EXPECT_NEAR(state.mu, gold1.mu, 1e-8);
    EXPECT_NEAR(state.k, gold1.k, 1e-4);
    EXPECT_NEAR(state.v, gold1.v, 1e-9);
    EXPECT_NEAR(state.s, gold1.s, 3.17);
    EXPECT_NEAR(state.h, gold1.h, 132);
    EXPECT_NEAR(state.w, gold1.w, 0.12);
}

TEST(InterpolatedFluidPropertiesTest, h_s)
{
    auto file_name = path(FPROPS_UNIT_TESTS_ROOT) / path("assets") / path("n2-smol.fprops.h5");
    InterpolatedFluidProperties fp;
    fp.load(file_name.string());

    double h = 288084;
    double s = 6083;
    EXPECT_THROW({ fp.h_s(h, s); }, std::runtime_error);
}

TEST(InterpolatedFluidPropertiesTest, non_existent_file)
{
    InterpolatedFluidProperties fp;
    EXPECT_THROW({ fp.load("non-existent-file"); }, std::runtime_error);
}

TEST(InterpolatedFluidPropertiesTest, empty_file)
{
    auto file_name = path(FPROPS_UNIT_TESTS_ROOT) / path("assets") / path("empty.fprops.h5");
    InterpolatedFluidProperties fp;
    fp.load(file_name);
    EXPECT_THROW({ fp.p_T(1e6, 280); }, std::runtime_error);
    EXPECT_THROW({ fp.rho_T(0.1, 280); }, std::runtime_error);
    EXPECT_THROW({ fp.rho_p(0.1, 1e6); }, std::runtime_error);
    EXPECT_THROW({ fp.v_u(1, 1e5); }, std::runtime_error);
    EXPECT_THROW({ fp.h_s(1e4, 1.1e4); }, std::runtime_error);
}

TEST(InterpolatedFluidPropertiesTest, grid_1_by_1_file)
{
    auto file_name = path(FPROPS_UNIT_TESTS_ROOT) / path("assets") / path("grid-1x1.fprops.h5");
    InterpolatedFluidProperties fp;
    EXPECT_THROW({ fp.load(file_name); }, std::runtime_error);
}

TEST(InterpolatedFluidPropertiesTest, wrong_units_file)
{
    auto file_name = path(FPROPS_UNIT_TESTS_ROOT) / path("assets") / path("wrong-units.fprops.h5");
    InterpolatedFluidProperties fp;
    EXPECT_THROW({ fp.load(file_name); }, std::runtime_error);
}
