cmake_minimum_required(VERSION 2.8)
project(samtools)

file(GLOB SOURCE_FILES *.c)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/samtools/htslib/)
add_library(samtools ${SOURCE_FILES})
