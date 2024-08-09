// SPDX-FileCopyrightText: 2024 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#pragma once

#include "fmt/format.h"
#include <exception>

namespace fprops {

class Exception : public std::exception {
public:
    template <typename... T>
    Exception(fmt::format_string<T...> format, T... args) :
        msg(fmt::format(format, std::forward<T>(args)...))
    {
    }

    /// Get the exception message
    [[nodiscard]] auto what() const noexcept -> const char * override;

private:
    /// Error message
    std::string msg;
};

} // namespace fprops
