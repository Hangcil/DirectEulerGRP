cmake_minimum_required(VERSION 3.5.0)
project(GRPForEulerEq LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /fp:precise")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /GL")
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /LTCG")
set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} /LTCG")



find_package(Eigen3 REQUIRED PATHS "C:\\libs\\Eigen3\\share\\eigen3\\cmake")
include_directories(${EIGEN3_INCLUDE_DIRS})
add_executable(GRPForEulerEq ExampleLeblanc.cpp RPSolver.cpp RPSolver.h GRPSolver.h GRPSolver.cpp FVM.h FVM.cpp)



