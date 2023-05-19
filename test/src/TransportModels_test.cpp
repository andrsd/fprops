#include "gtest/gtest.h"
#include "TransportModels.h"

using namespace fprops;

TEST(TransportTest, eta0_and_poly)
{
    Eta0AndPoly<double> eta0p({ 2, 3 }, { 3, 4 });
    EXPECT_DOUBLE_EQ(eta0p.value(12., 2.), 72.);
}

TEST(TransportTest, lennard_jones)
{
    LennardJones<double> lj(2.e-3, 3., 2., { 2, 3 });
    EXPECT_DOUBLE_EQ(lj.value(2.), 0.0060967411665096916);
}

TEST(TransportTest, modified_batshinski_hildebrand)
{
    ModifiedBatshinskiHildebrand<double> mbh({ 2, 3 }, { 3, 4 }, { 4, 5 }, { 4, 3 }, { 2, 3 });
    EXPECT_DOUBLE_EQ(mbh.value(2., 3.), 0.000097523945419603);
}
