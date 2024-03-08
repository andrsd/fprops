#include "gtest/gtest.h"
#include "ExceptionTestMacros.h"
#include "fprops/Numerics.h"

using namespace fprops;

TEST(NumericsTest, sqr)
{
    EXPECT_DOUBLE_EQ(sqr(2.), 4.);
}

TEST(NumericsTest, cb)
{
    EXPECT_DOUBLE_EQ(cb(2.), 8.);
}

TEST(NumericsTest, pow)
{
    EXPECT_DOUBLE_EQ(math::pow(2., 3), 8.);
    EXPECT_DOUBLE_EQ(math::pow(3., 3), 27.);
    EXPECT_DOUBLE_EQ(math::pow(5., 4), 625.);

    EXPECT_DOUBLE_EQ(math::pow(2., -1), 0.5);
}

TEST(NumericsTest, newton_root)
{
    auto f = [](double x) {
        return x * x - 4;
    };

    auto df = [](double x) {
        return 2 * x;
    };

    EXPECT_DOUBLE_EQ(newton::root(3., f, df), 2);
    EXPECT_DOUBLE_EQ(newton::root(-3., f, df), -2);
}

TEST(NumericsTest, newton_root_diverge)
{
    auto f = [](double x) {
        return x * x * x - 2 * x + 2;
    };

    auto df = [](double x) {
        return 3 * x * x - 2;
    };

    EXPECT_THROW_MSG(newton::root(0, f, df), "Newton's method failed to converge");
}
