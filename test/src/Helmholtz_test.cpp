#include "gmock/gmock.h"
#include "Helmholtz.h"

using namespace fprops;

namespace {

class MockHelmholtz : public Helmholtz {
public:
    MockHelmholtz() : Helmholtz(1., 1., 1., 1.) {}

    MOCK_METHOD(double, alpha, (double delta, double tau), (const));
    MOCK_METHOD(double, dalpha_ddelta, (double delta, double tau), (const));
    MOCK_METHOD(double, dalpha_dtau, (double delta, double tau), (const));
    MOCK_METHOD(double, d2alpha_ddelta2, (double delta, double tau), (const));
    MOCK_METHOD(double, d2alpha_dtau2, (double delta, double tau), (const));
    MOCK_METHOD(double, d2alpha_ddeltatau, (double delta, double tau), (const));

    MOCK_METHOD(double, mu_from_rho_T, (double rho, double T), (const));
    MOCK_METHOD(double, k_from_rho_T, (double rho, double T), (const));
};

} // namespace

TEST(HelmholtzTest, p_T_incorrect)
{
    MockHelmholtz fp;

    EXPECT_THROW(auto p = fp.p_T(1e5, -1), std::domain_error);
}

TEST(HelmholtzTest, v_u_incorrect)
{
    MockHelmholtz fp;

    EXPECT_THROW(auto p = fp.v_u(-1, 1), std::domain_error);
    EXPECT_THROW(auto p = fp.v_u(1, -1), std::domain_error);
}

TEST(HelmholtzTest, ideal_gas_lead)
{
    class Test : public MockHelmholtz {
    public:
        Test() : MockHelmholtz(), a(10., 100.)
        {
            double e = std::exp(1.);
            EXPECT_DOUBLE_EQ(a.alpha(e, 1), 111.);
            EXPECT_DOUBLE_EQ(a.ddelta(10, 0), 0.1);
            EXPECT_DOUBLE_EQ(a.dtau(0, 0), 100.);
            EXPECT_DOUBLE_EQ(a.d2delta(2, 0), -0.25);
        }

        Helmholtz::IdealGasLead<double> a;
    };

    Test t;
}

TEST(HelmholtzTest, ideal_gas_log_tau)
{
    class Test : public MockHelmholtz {
    public:
        Test() : MockHelmholtz(), a(11.)
        {
            double e = std::exp(1.);
            EXPECT_DOUBLE_EQ(a.alpha(1., e), 11.);
            EXPECT_DOUBLE_EQ(a.dtau(0, 2), 5.5);
            EXPECT_DOUBLE_EQ(a.d2tau(0, 2), -2.75);
        }

        Helmholtz::IdealGasLogTau<double> a;
    };

    Test t;
}

TEST(HelmholtzTest, ideal_gas_power)
{
    class Test : public MockHelmholtz {
    public:
        Test() : MockHelmholtz(), a({ 3, 4 }, { 2, 3 })
        {
            EXPECT_DOUBLE_EQ(a.alpha(0., 2), 44.);
            EXPECT_DOUBLE_EQ(a.dtau(0, 2), 60);
            EXPECT_DOUBLE_EQ(a.d2tau(0, 2), 54);
        }

        Helmholtz::IdealGasPower<double> a;
    };

    Test t;
}

TEST(HelmholtzTest, enthalpy_entropy_offset)
{
    class Test : public MockHelmholtz {
    public:
        Test() : MockHelmholtz(), a(11, 12)
        {
            EXPECT_DOUBLE_EQ(a.alpha(0., 2), 35.);
            EXPECT_DOUBLE_EQ(a.dtau(0, 2), 12.);
        }

        Helmholtz::IdealEnthalpyEntropyOffset<double> a;
    };

    Test t;
}

TEST(HelmholtzTest, ideal_gas_planck_einstein_generalized)
{
    class Test : public MockHelmholtz {
    public:
        Test() : MockHelmholtz(), a({ 2, 3 }, { 3, 2 }, { 3, 4 }, { 5, 6 })
        {
            EXPECT_DOUBLE_EQ(a.alpha(0., 1.), 20.9121719702137190);
            EXPECT_DOUBLE_EQ(a.dtau(0., 1.), 11.3294239412548978);
            EXPECT_DOUBLE_EQ(a.d2tau(0., 1.), 1.41785827094037168);
        }

        Helmholtz::IdealPlanckEinsteinGeneralized<double> a;
    };

    Test t;
}

TEST(HelmholtzTest, ideal_gas_planck_einstein_function_t)
{
    class Test : public MockHelmholtz {
    public:
        Test() : MockHelmholtz(), a(2, { 2, 3 }, { 3, 4 })
        {
            EXPECT_DOUBLE_EQ(a.alpha(0., 1.), -0.941205291457485162);
            EXPECT_DOUBLE_EQ(a.dtau(0., 1.), 1.800756606864598643);
            EXPECT_DOUBLE_EQ(a.d2tau(0., 1.), -3.835882116252504991);
        }

        Helmholtz::IdealGasPlanckEinsteinFunctionT<double> a;
    };

    Test t;
}

TEST(HelmholtzTest, residual_power)
{
    class Test : public MockHelmholtz {
    public:
        Test() : MockHelmholtz(), a({ 2, 3 }, { 4, 5 }, { 6, 7 })
        {
            EXPECT_DOUBLE_EQ(a.alpha(2., 3.), 233280);
            EXPECT_DOUBLE_EQ(a.ddelta(2., 3.), 571536);
            EXPECT_DOUBLE_EQ(a.dtau(2., 3.), 536544);
            EXPECT_DOUBLE_EQ(a.d2delta(2., 3.), 1119744);
            EXPECT_DOUBLE_EQ(a.d2tau(2., 3.), 1057536);
            EXPECT_DOUBLE_EQ(a.d2deltatau(2., 3.), 1318032);
        }

        Helmholtz::ResidualPower<double> a;
    };

    Test t;
}

TEST(HelmholtzTest, residual_power_exp)
{
    class Test : public MockHelmholtz {
    public:
        Test() : MockHelmholtz(), a({ 2, 3 }, { 3, 4 }, { 4, 5 }, { 2, 3 })
        {
            EXPECT_DOUBLE_EQ(a.alpha(2., 3.), 27.649904091654395748);
            EXPECT_DOUBLE_EQ(a.ddelta(2., 3.), -98.471030918047725031);
            EXPECT_DOUBLE_EQ(a.dtau(2., 3.), 38.170817486157493694);
            EXPECT_DOUBLE_EQ(a.d2delta(2., 3.), 423.496477990674375469);
            EXPECT_DOUBLE_EQ(a.d2tau(2., 3.), 40.344615314965770409);
            EXPECT_DOUBLE_EQ(a.d2deltatau(2., 3.), -144.337494863579960335);
        }

        Helmholtz::ResidualPowerExp<double, unsigned int> a;
    };

    Test t;
}

TEST(HelmholtzTest, residual_gaussian)
{
    class Test : public MockHelmholtz {
    public:
        Test() :
            MockHelmholtz(),
            a({ 2, 3 }, { 4, 5 }, { 6, 7 }, { 2, 3 }, { 4, 2 }, { 3, 2 }, { 2, 3 })
        {
            EXPECT_DOUBLE_EQ(a.alpha(2., 3.), 209952.389617276034850740);
            EXPECT_DOUBLE_EQ(a.ddelta(2., 3.), 524883.896172760348507404);
            EXPECT_DOUBLE_EQ(a.dtau(2., 3.), 489886.441530895860597038);
            EXPECT_DOUBLE_EQ(a.d2delta(2., 3.), -209914.986358776689179657);
            EXPECT_DOUBLE_EQ(a.d2tau(2., 3.), 139971.636427909658606910);
            EXPECT_DOUBLE_EQ(a.d2deltatau(2., 3.), 1.224704415308958605e6);
        }

        Helmholtz::ResidualGaussian<double> a;
    };

    Test t;
}
