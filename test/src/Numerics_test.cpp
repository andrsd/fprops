#include "gtest/gtest.h"
#include "Numerics.h"

using namespace fprops;

TEST(NumericsTest, sqr)
{
    EXPECT_DOUBLE_EQ(sqr(2.), 4.);
}

TEST(NumericsTest, cb)
{
    EXPECT_DOUBLE_EQ(cb(2.), 8.);
}
