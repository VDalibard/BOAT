cmake_minimum_required (VERSION 3.5)
project(bo_naive)

set(CMAKE_CXX_COMPILER g++-8)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "-Wall -fexceptions -O3 -fopenmp")
# set(CMAKE_VERBOSE_MAKEFILE on)

set(BOAT_DIR "./../..")

include_directories(SYSTEM $ENV{EIGEN_DIR})

include_directories("${BOAT_DIR}/include/")
include_directories("${BOAT_DIR}/build/")

find_library(NLOPT_LIBRARY nlopt)
find_library(BOAT_LIBRARY boat PATHS ./../../build/)

add_executable(bo_naive_opti bo_naive_opti.cpp)

target_link_libraries(bo_naive_opti pthread blas nlopt boost_serialization boost_system boost_coroutine ${BOAT_LIBRARY})