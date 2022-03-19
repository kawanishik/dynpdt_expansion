cmake_minimum_required(VERSION 3.17)

project(check CXX)

add_executable(check check.cpp)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_library(check INTERFACE poplar)

