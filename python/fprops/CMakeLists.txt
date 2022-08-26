cmake_minimum_required(VERSION 3.4...3.18)
project(pyfprops LANGUAGES CXX)

find_package(pybind11 REQUIRED)

pybind11_add_module(pyfprops fprops.cpp)

target_include_directories(
    ${PROJECT_NAME}
    PRIVATE
        ${CMAKE_BINARY_DIR}
        ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(
    ${PROJECT_NAME}
    PUBLIC
        fprops
)
