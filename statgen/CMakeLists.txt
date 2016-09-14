cmake_minimum_required(VERSION 2.8)
project(statgen)

file(GLOB CXX_SOURCE_FILES *.cpp)
file(GLOB C_SOURCE_FILES *.c)

add_library(statgen ${CXX_SOURCE_FILES} ${C_SOURCE_FILES})
