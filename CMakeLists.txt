cmake_minimum_required(VERSION 3.19)

project(
    fprops
    VERSION 8.0.0
    LANGUAGES C CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

include(FetchContent)
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)
include(CodeCoverage)
include(Sanitization)

option(FPROPS_WITH_PYTHON "Build python wrapper" YES)
option(FPROPS_BUILD_TESTS "Build tests" NO)

set(HDF5_PREFER_PARALLEL YES)
find_package(HDF5 1.8 COMPONENTS C HL REQUIRED)
if (HDF5_IS_PARALLEL)
    find_package(MPI REQUIRED)
endif()
find_package(Eigen3 3.3.4 REQUIRED)
find_package(fmt 11 REQUIRED)
if (FPROPS_WITH_PYTHON)
    add_subdirectory(python)
endif()

if (FPROPS_BUILD_TESTS)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        v1.15.2
    )
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    mark_as_advanced(FORCE
        BUILD_GMOCK
        INSTALL_GTEST
        GTEST_HAS_ABSL
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
endif()

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
