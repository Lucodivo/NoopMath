cmake_minimum_required(VERSION 3.28)
project("noop_math")

set(CMAKE_CXX_STANDARD 17)

add_library(noop_math STATIC noop_math.cpp)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tests ${CMAKE_CURRENT_SOURCE_DIR/tests/build})