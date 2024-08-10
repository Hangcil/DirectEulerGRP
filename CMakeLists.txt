cmake_minimum_required(VERSION 3.5.0)
project(GRPForEulerEq LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED PATHS "C:\\libs\\Eigen3\\share\\eigen3\\cmake")
include_directories(${EIGEN3_INCLUDE_DIRS} src)
add_executable(GRPForEulerEq Example9_1.cpp src/FVM.h src/FVM.cpp src/solvers.h src/solvers.cpp)


