project(fprops-test)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(src)

add_test(
    NAME fprops-test
    COMMAND ${PROJECT_NAME}
)
