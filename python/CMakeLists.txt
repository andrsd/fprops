find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
find_program(PYTHON "python" REQUIRED)

configure_file(setup.py.in setup.py)
configure_file(src/fprops/__init__.py src/fprops/__init__.py)
configure_file(src/fprops/_version.py.in src/fprops/_version.py)
install(CODE "execute_process(COMMAND ${PYTHON} -m pip install ${CMAKE_CURRENT_BINARY_DIR})")

add_subdirectory(src)

if (FPROPS_BUILD_TESTS)
    find_program(PYTEST "pytest" REQUIRED)
    enable_testing()
    add_test(
        NAME python-tests
        COMMAND ${PYTEST} ${CMAKE_CURRENT_SOURCE_DIR}
    )
    set_tests_properties(python-tests PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}/src:$ENV{PYTHONPATH}")
endif()
