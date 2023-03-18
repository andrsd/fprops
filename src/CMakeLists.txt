add_library(
    ${PROJECT_NAME}
    SHARED
        FluidProperties.cpp
        Helmholtz.cpp
        IdealGas.cpp
        Nitrogen.cpp
        Numerics.cpp
        SinglePhaseFluidProperties.cpp
)
if(CMAKE_PROJECT_NAME STREQUAL "fprops")
    target_code_coverage(${PROJECT_NAME})
    target_sanitization(${PROJECT_NAME})
endif()

set_target_properties(
    ${PROJECT_NAME}
    PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}
)

target_include_directories(
    ${PROJECT_NAME}
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/fprops>
    PRIVATE
        ${CMAKE_BINARY_DIR}
)

target_link_libraries(
    ${PROJECT_NAME}
    PUBLIC
        fmt::fmt
)

install(
    TARGETS fprops
    EXPORT fpropsTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
    DIRECTORY ${CMAKE_SOURCE_DIR}/include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/fprops
    FILES_MATCHING PATTERN "*.h"
)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    fpropsConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(
    EXPORT fpropsTargets
    FILE fpropsTargets.cmake
    NAMESPACE fprops::
    DESTINATION lib/cmake/fprops
)

install(
    FILES
        "${PROJECT_SOURCE_DIR}/cmake/fpropsConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/fpropsConfigVersion.cmake"
    DESTINATION
        lib/cmake/${PROJECT_NAME}
)
