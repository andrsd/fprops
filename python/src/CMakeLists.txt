cmake_minimum_required(VERSION 3.4...3.18)
project(pyfprops LANGUAGES C CXX)

find_package(pybind11 2.9 REQUIRED)

pybind11_add_module(${PROJECT_NAME} fprops.cpp)

target_include_directories(
    ${PROJECT_NAME}
    PRIVATE
        ${CMAKE_BINARY_DIR}
        ${CMAKE_SOURCE_DIR}/include
        ${EIGEN3_INCLUDE_DIR}
)

target_link_libraries(
    ${PROJECT_NAME}
    PUBLIC
        fprops
)

install(
    TARGETS ${PROJECT_NAME}
    COMPONENT python
    LIBRARY DESTINATION "lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages/pyfprops"
)
