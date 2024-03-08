#include "fprops/Exception.h"

namespace fprops {

auto
Exception::what() const noexcept -> const char *
{
    return this->msg.c_str();
}

} // namespace fprops
