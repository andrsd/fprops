#pragma once

#include "fprops/Exception.h"

/// Test that the `cmd` will throw `fprops::Exception` with message `msg`
#define EXPECT_THROW_MSG(cmd, msg)   \
    try {                            \
        cmd;                         \
        FAIL();                      \
    }                                \
    catch (Exception & e) {          \
        EXPECT_STREQ(e.what(), msg); \
    }                                \
    catch (...) {                    \
        FAIL();                      \
    }

#define EXPECT_THAT_THROW_MSG(cmd, matcher) \
    try {                                   \
        cmd;                                \
        FAIL();                             \
    }                                       \
    catch (Exception & e) {                 \
        EXPECT_THAT(e.what(), matcher);     \
    }                                       \
    catch (...) {                           \
        FAIL();                             \
    }
