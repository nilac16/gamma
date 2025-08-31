cmake_minimum_required(VERSION 3.15)

set(CMAKE_USER_MAKE_RULES_OVERRIDE "cmake/BuildFlags.cmake")

set(CMAKE_C_STANDARD 99)

project(gamma C CXX)

find_package(DCMTK REQUIRED)

add_subdirectory(gamma)
