cmake_minimum_required(VERSION 2.8)
project(samtools)

file(GLOB SOURCE_FILES *.c)
include_directories(${samtools_SOURCE_DIR}/htslib/)
add_library(samtools ${SOURCE_FILES})
