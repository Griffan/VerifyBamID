cmake_minimum_required(VERSION 2.8)
project(libVcf)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "/Users/fanzhang/BioLib")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

file(GLOB SOURCE_FILES
    "*.h"
    "*.cpp"
)
include_directories(${LIBSTATGEN_INCLUDE_DIR})
add_library(Vcf ${SOURCE_FILES})
#set_target_properties(Vcf PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/libVcf")
target_link_libraries(Vcf statgen)
#add_dependencies(Vcf statgen )

