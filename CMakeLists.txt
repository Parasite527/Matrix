cmake_minimum_required(VERSION 3.20)
project(Matrix)

set(CMAKE_CXX_STANDARD 17)

set(GCC_COVERAGE_COMPILE_FLAGS "-Wall -Wextra -Wconversion -fsanitize==address")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )

add_executable(Matrix main.cpp)
