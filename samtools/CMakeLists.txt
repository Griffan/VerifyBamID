cmake_minimum_required(VERSION 2.8)
project(samtools)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3  -std=c++11")

file(GLOB SOURCE_FILES *.c)

add_library(samtools ${SOURCE_FILES})
