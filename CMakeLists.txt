cmake_minimum_required(VERSION 3.15)

set(CMAKE_USER_MAKE_RULES_OVERRIDE "cmake/BuildFlags.cmake")

set(CMAKE_C_STANDARD 11)

project(gamma C)

add_subdirectory(gamma)
