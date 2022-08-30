find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
find_program(PYTHON "python")

configure_file(setup.py.in setup.py)
install(CODE "execute_process(COMMAND ${PYTHON} -m pip install ${CMAKE_CURRENT_BINARY_DIR})")

add_subdirectory(src)
