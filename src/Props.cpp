#include "fprops/State.h"

namespace fprops {

State::State() :
    u(0.),
    v(0.),
    rho(0.),
    p(0.),
    T(0.),
    mu(0.),
    cp(0.),
    cv(0.),
    s(0.),
    k(0.),
    h(0.),
    w(0.)
{
}

State::State(double u,
             double v,
             double rho,
             double p,
             double T,
             double mu,
             double cp,
             double cv,
             double s,
             double k,
             double h,
             double w) :
    u(u),
    v(v),
    rho(rho),
    p(p),
    T(T),
    mu(mu),
    cp(cp),
    cv(cv),
    s(s),
    k(k),
    h(h),
    w(w)
{
}

} // namespace fprops
