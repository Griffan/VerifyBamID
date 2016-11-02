cmake_minimum_required(VERSION 2.8)
project(VerifyBamID)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3  -std=c++0x -Wno-unused-variable -Wno-unused-result -Wno-format-security")
add_definitions(-D__STDC_LIMIT_MACROS)

include(ExternalProject)
ExternalProject_Add(HTSLIB
        PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
        STAMP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/build/stamp
        GIT_REPOSITORY https://github.com/samtools/htslib.git
        GIT_TAG "master"

        UPDATE_COMMAND ""

        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/samtools/htslib
        CONFIGURE_COMMAND autoheader && autoconf && ${CMAKE_CURRENT_SOURCE_DIR}/samtools/htslib/configure --prefix=${CMAKE_CURRENT_SOURCE_DIR}/samtools/htslib
        BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/samtools/htslib
        BUILD_COMMAND make && make test
)

add_subdirectory(statgen)
add_subdirectory(samtools)


find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(SOURCE_FILES main.cpp ContaminationEstimator.cpp MathGenMin.cpp MathGold.cpp Random.cpp SimplePileupViewer.cpp params.cpp )#SVDcalculator.cpp SVDcalculator.h)
include_directories(statgen ${CMAKE_CURRENT_SOURCE_DIR}/samtools/htslib/include/ samtools)
add_executable(VerifyBamID ${SOURCE_FILES})

find_library(ZLIB z HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
#find_library(HITSLIB hts HINTS  ${CMAKE_CURRENT_SOURCE_DIR}/samtools/htslib/ /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)

target_link_libraries(VerifyBamID statgen samtools ${ZLIB})
add_dependencies(VerifyBamID HTSLIB)
enable_testing()
add_test( NAME myTest
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
          COMMAND VerifyBamID --UDPath ./resource/hapmap_3.3.b37.dat.UD --BamFile ./resource/test.bam --BedPath ./resource/choose.bed --MeanPath ./resource/hapmap_3.3.b37.dat.mu --Reference ./resource/chr20.fa.gz )