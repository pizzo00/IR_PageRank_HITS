cmake_minimum_required(VERSION 3.22)
project(IR_Projects)

set(CMAKE_CXX_STANDARD 20)

#set(CMAKE_CXX_FLAGS "-O3")

add_executable(IR_Projects main.cpp
        Sorter.cpp
        Sorter.h)
