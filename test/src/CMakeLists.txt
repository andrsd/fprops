add_executable(
    ${PROJECT_NAME}
        Air_test.cpp
        CarbonDioxide_test.cpp
        Helium_test.cpp
        Helmholtz_test.cpp
        IdealGas_test.cpp
        Nitrogen_test.cpp
        Numerics_test.cpp
        InterpolatedFluidProperties_test.cpp
        TransportModels_test.cpp
        Utils_test.cpp
        main.cpp
)
target_code_coverage(${PROJECT_NAME})
target_sanitization(${PROJECT_NAME})
target_compile_definitions(${PROJECT_NAME} PRIVATE H5PP_USE_EIGEN3)
add_compile_definitions("FPROPS_UNIT_TESTS_ROOT=\"${FPROPS_UNIT_TESTS_ROOT}\"")

target_include_directories(
    ${PROJECT_NAME}
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
        ${CMAKE_BINARY_DIR}
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/contrib
    PRIVATE
        ${EIGEN3_INCLUDE_DIR}
)

target_link_libraries(
    ${PROJECT_NAME}
    PUBLIC
        libfprops
        gtest_main
        gmock_main
)
