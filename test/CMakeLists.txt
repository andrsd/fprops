project(fprops-test)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(src)

add_test(
    NAME fprops-test
    COMMAND ${PROJECT_NAME}
)
set_tests_properties(${PROJECT_NAME}
    PROPERTIES
        LABELS unit-tests
)

if(FPROPS_CODE_COVERAGE)
    set_tests_properties(
        ${PROJECT_NAME}
        PROPERTIES
            ENVIRONMENT LLVM_PROFILE_FILE=${PROJECT_NAME}.profraw
    )
endif()
