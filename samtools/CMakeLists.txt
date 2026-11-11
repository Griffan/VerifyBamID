cmake_minimum_required(VERSION 3.10)
project(samtools)

file(GLOB SOURCE_FILES *.c)

add_library(samtools ${SOURCE_FILES})
include_directories(samtools)
target_link_htslib(samtools)
