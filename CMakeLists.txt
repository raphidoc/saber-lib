cmake_minimum_required(VERSION 3.10)
project(saber-lib C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Header and source files
include_directories(include)
file(GLOB SRC_FILES "src/*.c")

# Build static library
add_library(saber STATIC ${SRC_FILES})

# Build test executable
add_executable(saber_test test/main.c)
target_link_libraries(saber_test saber m)