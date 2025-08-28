cmake_minimum_required(VERSION 3.15)

set(CMAKE_USER_MAKE_RULES_OVERRIDE "cmake/BuildFlags.cmake")

project(gamma C CXX)

find_package(DCMTK REQUIRED)

add_subdirectory(gamma)
