cmake_minimum_required(VERSION 3.5.0)
project(GRPForEulerEq LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED PATHS "C:\\libs\\Eigen3\\share\\eigen3\\cmake")
include_directories(${EIGEN3_INCLUDE_DIRS} src)
add_executable(GRPForEulerEq src/basicFuns.h src/basicFuns.cpp src/wrappers.h src/wrappers.cpp src/solvers.h src/solvers.cpp "TestAccuracy.cpp")


