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

TEST(Helmholtz, p_T_incorrect)
{
    MockHelmholtz fp;

    EXPECT_THROW(auto p = fp.p_T(1e5, -1), std::domain_error);
}

TEST(Helmholtz, v_u_incorrect)
{
    MockHelmholtz fp;

    EXPECT_THROW(auto p = fp.v_u(-1, 1), std::domain_error);
    EXPECT_THROW(auto p = fp.v_u(1, -1), std::domain_error);
}
