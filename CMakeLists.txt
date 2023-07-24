cmake_minimum_required(VERSION VERSION 3.16)

project(
    fprops
    VERSION 1.0
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

include(FetchContent)
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)
include(${PROJECT_SOURCE_DIR}/cmake/CodeCoverage.cmake)
include(${PROJECT_SOURCE_DIR}/cmake/Sanitization.cmake)

option(FPROPS_WITH_PYTHON "Build python wrapper" YES)
option(FPROPS_BUILD_TESTS "Build tests" NO)

find_package(fmt REQUIRED)
if (FPROPS_WITH_PYTHON)
    add_subdirectory(python)
endif()

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        release-1.11.0
)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
mark_as_advanced(FORCE
    BUILD_GMOCK
    INSTALL_GTEST
    FETCHCONTENT_SOURCE_DIR_GOOGLETEST
    FETCHCONTENT_UPDATES_DISCONNECTED_GOOGLETEST
)

FetchContent_MakeAvailable(googletest)
mark_as_advanced(FORCE
    FETCHCONTENT_BASE_DIR
    FETCHCONTENT_FULLY_DISCONNECTED
    FETCHCONTENT_QUIET
    FETCHCONTENT_UPDATES_DISCONNECTED
)


add_subdirectory(src)

# Tests

if (FPROPS_BUILD_TESTS)
    set(FPROPS_UNIT_TESTS_ROOT ${CMAKE_SOURCE_DIR}/test CACHE PATH "" FORCE)
    mark_as_advanced(FORCE
        FPROPS_UNIT_TESTS_ROOT
    )

    enable_testing()
    add_subdirectory(test)
endif()

add_subdirectory(docs)
