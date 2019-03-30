cmake_minimum_required (VERSION 3.0)
project (CS888Project)

include_directories (include)
file (GLOB SOURCES "src/*.cpp")

add_executable(CS888Project ${SOURCES})
