#include "gtest/gtest.h"
#include "fprops/transport_models.h"

using namespace fprops;

TEST(TransportTest, polynomial)
{
    Polynomial<double> poly({ 2, 3 }, { 3, 4 }, { 5, 6 });
    EXPECT_DOUBLE_EQ(poly.value(12., 2.), 4091904.);
}

TEST(TransportTest, eta0_and_poly)
{
    Eta0AndPoly<double> eta0p({ 2, 3 }, { 3, 4 });
    EXPECT_DOUBLE_EQ(eta0p.value(12., 2.), 72.);
}

TEST(TransportTest, lennard_jones)
{
    LennardJones<double> lj(0.0266958, 2.e-3, 3., 2., { 2, 3 });
    EXPECT_DOUBLE_EQ(lj.value(2.), 0.0060967411665096916);
}

TEST(TransportTest, modified_batshinski_hildebrand)
{
    ModifiedBatshinskiHildebrand<double> mbh({ 2, 3 }, { 3, 4 }, { 4, 5 }, { 4, 3 }, { 2, 3 });
    EXPECT_DOUBLE_EQ(mbh.value(2., 3.), 0.000097523945419603);
}

TEST(TransportTest, powers_of_T)
{
    PowersOfTemperature<double> pot({ 3. }, { 5. });
    EXPECT_DOUBLE_EQ(pot.value(2), 96.);
}

TEST(TransportTest, powers_of_T_reduced)
{
    PowersOfTreduced<double> potr(300, { 3. }, { 5. });
    EXPECT_DOUBLE_EQ(potr.value(600), 96.);
}

TEST(TransportTest, rainwater_friend)
{
    RainwaterFriend<double> rf(10, 3e-5, { 4, 5 }, { 2, 3 });
    EXPECT_DOUBLE_EQ(rf.value(25), 1676789965434.375);
}

TEST(TransportTest, quadratic_general_ft)
{
    QuadraticGeneralFrictionTheory<double> qgft({ 2, 3, 4 },
                                                { 3, 4, 5 },
                                                { 5, 6, 7 },
                                                { 1, 3, 5 },
                                                { 5, 7, 9 },
                                                { 11, 17, 3 });
    EXPECT_DOUBLE_EQ(qgft.value(1e-4, 1.5, 2., 5.4), 0.0019802334193187326);
}
