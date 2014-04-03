cmake_minimum_required(VERSION 2.8.0) 
project(mpdata C CXX Fortran)

enable_testing()
add_subdirectory(cpp)
add_subdirectory(for)
add_subdirectory(tests)
