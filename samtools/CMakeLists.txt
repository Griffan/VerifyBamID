cmake_minimum_required(VERSION 2.8)
project(samtools)

file(GLOB SOURCE_FILES *.c)

add_library(samtools ${SOURCE_FILES})
include_directories(samtools ${HTS_INCLUDE_DIRS})
target_link_libraries(samtools ${HTS_LIBRARIES})
