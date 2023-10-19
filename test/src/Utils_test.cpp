#include "gmock/gmock.h"
#include "Utils.h"

using namespace fprops;

TEST(UtilsTest, interval_index)
{
    std::vector<double> range({ 1., 2.5, 7., 9.6, 12. });
    EXPECT_EQ(utils::interval_index(range, 1.), 0);
    EXPECT_EQ(utils::interval_index(range, 2.), 0);
    EXPECT_EQ(utils::interval_index(range, 2.5), 1);
    EXPECT_EQ(utils::interval_index(range, 5.), 1);
    EXPECT_EQ(utils::interval_index(range, 7.25), 2);
    EXPECT_EQ(utils::interval_index(range, 10), 3);
}

TEST(UtilsTest, normalize_interval_location)
{
    EXPECT_DOUBLE_EQ(utils::normalize_interval_location(1., 3., 1.), 0.);
    EXPECT_DOUBLE_EQ(utils::normalize_interval_location(1., 3., 2.), 0.5);
    EXPECT_DOUBLE_EQ(utils::normalize_interval_location(1., 3., 3.), 1.);

    // test on interval of zero length
    EXPECT_DOUBLE_EQ(utils::normalize_interval_location(1., 1., 1.), 0.);
}
