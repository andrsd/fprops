find_package(Python3 REQUIRED COMPONENTS Interpreter Development)

configure_file(pyproject.toml.in pyproject.toml)
install(CODE "execute_process(COMMAND ${Python3_EXECUTABLE} -m pip install ${CMAKE_CURRENT_BINARY_DIR})")

add_subdirectory(src)

if (FPROPS_BUILD_TESTS)
    find_program(PYTEST "pytest" REQUIRED)
    enable_testing()
    add_test(
        NAME python-tests
        COMMAND ${PYTEST} ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()
