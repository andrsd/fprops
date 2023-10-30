project(fprops LANGUAGES C CXX)

find_package(pybind11 2.9 REQUIRED)

pybind11_add_module(fprops fprops.cpp)

target_include_directories(
    fprops
    PRIVATE
        ${CMAKE_BINARY_DIR}
        ${CMAKE_SOURCE_DIR}/include
        ${EIGEN3_INCLUDE_DIR}
)

target_link_libraries(
    fprops
    PUBLIC
        libfprops
)

install(
    TARGETS fprops
    COMPONENT python
    LIBRARY DESTINATION "lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages/fprops"
)
