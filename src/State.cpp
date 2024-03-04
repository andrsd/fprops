#include "fprops/State.h"
#include "fmt/printf.h"

namespace fprops {

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
