cmake_minimum_required(VERSION 3.2)

project(NS_Solver_RK4)
add_compile_options(-std=c++20 -Wall -Wextra)

set(SOURCE_FILES main.cpp mesh.hpp)
add_executable(exec ${SOURCE_FILES})
