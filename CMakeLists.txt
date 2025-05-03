cmake_minimum_required(VERSION 3.10)
project(saber-lib
        VERSION 0.1.1
        LANGUAGES C
)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)

# Public headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# All source files
file(GLOB_RECURSE SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c")

# Build the library (static + shared if BUILD_SHARED_LIBS ON)
add_library(saber ${SRC_FILES})
target_link_libraries(saber PUBLIC m)

# Tell CMake which headers to install
set_target_properties(saber PROPERTIES
        PUBLIC_HEADER "include/saber.h"
        OUTPUT_NAME saber
)

# Generate the pkg-config file
configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/saber-lib.pc.in
        ${CMAKE_CURRENT_BINARY_DIR}/saber-lib.pc
        @ONLY
)

# Install rules
install(TARGETS saber
        EXPORT saberTargets
        ARCHIVE DESTINATION lib        # .a
        LIBRARY DESTINATION lib        # .so/.dylib
        PUBLIC_HEADER DESTINATION include/saber
)

install(EXPORT saberTargets
        NAMESPACE saber::
        DESTINATION lib/cmake/saber
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/saber-lib.pc
        DESTINATION lib/pkgconfig
)
