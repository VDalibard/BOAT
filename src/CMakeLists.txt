cmake_minimum_required(VERSION 2.8)
project(BOAT)

set(CMAKE_CXX_COMPILER g++-8)
set(CMAKE_CXX_FLAGS "-Wall -fexceptions -std=c++14 -O3 -fopenmp")
include_directories(SYSTEM $ENV{EIGEN_DIR})

add_subdirectory(optimization)
add_subdirectory(utilities)
add_subdirectory(probabilistic)
add_subdirectory(parameters)

add_library(boat $<TARGET_OBJECTS:optimization> $<TARGET_OBJECTS:gp> $<TARGET_OBJECTS:probabilistic> $<TARGET_OBJECTS:parameters>)

