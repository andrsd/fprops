project(fprops LANGUAGES C CXX)

find_package(pybind11 2.9 REQUIRED)

pybind11_add_module(fprops fprops.cpp)

target_include_directories(
    fprops
    PRIVATE
        ${CMAKE_BINARY_DIR}
        ${CMAKE_SOURCE_DIR}/include
        ${EIGEN3_INCLUDE_DIR}
        ${PROJECT_BINARY_DIR}
)

target_link_libraries(
    fprops
    PUBLIC
        libfprops
)

configure_file(version.h.in version.h)

set(PYTHON_SITE lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages)

install(
    TARGETS fprops
    COMPONENT python
    LIBRARY DESTINATION ${PYTHON_SITE}/fprops
)

install(
    FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/fprops/__init__.py
    DESTINATION ${PYTHON_SITE}/fprops
    COMPONENT python
)
