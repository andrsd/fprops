project(
    fprops
        VERSION 0.1
        LANGUAGES CXX
)

add_library(${PROJECT_NAME}
    SHARED
        BrentsMethod.cpp
        FluidProperties.cpp
        Helmholtz.cpp
        IdealGas.cpp
        Nitrogen.cpp
        SinglePhaseFluidProperties.cpp
)

set_target_properties(
    ${PROJECT_NAME}
    PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}
)

target_include_directories(
    ${PROJECT_NAME}
    PRIVATE
        ${CMAKE_BINARY_DIR}
        ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(
    ${PROJECT_NAME}
    PUBLIC
        fmt::fmt
)
