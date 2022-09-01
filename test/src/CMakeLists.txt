add_executable(${PROJECT_NAME}
    IdealGas_test.cpp
    Nitrogen_test.cpp
    Numerics_test.cpp
    main.cpp
)

target_include_directories(
    ${PROJECT_NAME}
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
        ${CMAKE_BINARY_DIR}
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/contrib
)

target_link_libraries(
    ${PROJECT_NAME}
    PUBLIC
        fprops
        gtest_main
        gmock_main
)

if (DEFINED ENV{CONDA_PREFIX})
    set_target_properties(${PROJECT_NAME} PROPERTIES BUILD_RPATH $ENV{CONDA_PREFIX}/lib)
endif()
