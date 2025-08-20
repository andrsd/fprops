// SPDX-FileCopyrightText: 2022 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#include "fprops/state.h"
#include "fmt/format.h"

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

auto
State::to_string() const -> std::string
{
    std::string str;
    str += fmt::format("rho = {} kg/m^3\n", this->rho);
    str += fmt::format("p = {} Pa\n", this->p);
    str += fmt::format("T = {} K\n", this->T);
    str += fmt::format("e = {} J/kg\n", this->u);
    str += fmt::format("v = {} m^3/kg\n", this->v);
    str += fmt::format("cp = {} J/(kg-K)\n", this->cp);
    str += fmt::format("cv = {} J/(kg-K)\n", this->cv);
    str += fmt::format("s = {} J/(kg-K)\n", this->s);
    str += fmt::format("h = {} J/kg\n", this->h);
    str += fmt::format("c = {} m/s\n", this->w);
    str += fmt::format("mu = {} Pa-s\n", this->mu);
    str += fmt::format("k = {} W/(m-K)\n", this->k);
    return str;
}

} // namespace fprops
