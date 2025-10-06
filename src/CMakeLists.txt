file(GLOB SRCS CONFIGURE_DEPENDS *.cpp)

add_library(libfprops SHARED ${SRCS})

target_compile_options(libfprops
    PRIVATE
        $<$<CONFIG:Debug>:-Wall -Werror>
)

if(NOT MSVC)
    set_target_properties(libfprops PROPERTIES OUTPUT_NAME fprops)
endif()
if(CMAKE_PROJECT_NAME STREQUAL "fprops")
    target_code_coverage(${PROJECT_NAME})
    target_sanitization(${PROJECT_NAME})
endif()

set_target_properties(
    libfprops
    PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}
)

target_include_directories(
    libfprops
    PUBLIC
        $<INSTALL_INTERFACE:include>
    PRIVATE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>/include
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>/include
        ${CMAKE_SOURCE_DIR}/contrib
        ${HDF5_INCLUDE_DIR}
        ${EIGEN3_INCLUDE_DIR}
)
if (HDF5_IS_PARALLEL)
    target_include_directories(
        libfprops
        PRIVATE
            ${MPI_C_INCLUDE_DIRS}
    )
endif()

target_link_libraries(
    libfprops
    PUBLIC
        fmt::fmt
    PRIVATE
        ${HDF5_LIBRARIES}
)

configure_package_config_file(
    ${PROJECT_SOURCE_DIR}/cmake/fprops-config.cmake.in fprops-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fprops
    NO_SET_AND_CHECK_MACRO
)
write_basic_package_version_file(fprops-config-version.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(
    TARGETS libfprops
    EXPORT fpropsTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
    DIRECTORY ${CMAKE_SOURCE_DIR}/include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.h"
)

install(
    EXPORT fpropsTargets
    FILE fprops-targets.cmake
    NAMESPACE fprops::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fprops
)

install(
    FILES
        "${CMAKE_CURRENT_BINARY_DIR}/fprops-config.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/fprops-config-version.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fprops
)
