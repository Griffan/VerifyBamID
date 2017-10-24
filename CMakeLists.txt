cmake_minimum_required(VERSION 2.8.4)
project(VerifyBamID)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3 -std=c++0x -Wno-unused-variable -Wno-unused-result -Wno-format-security")
add_definitions(-D__STDC_LIMIT_MACROS)
add_definitions(-D__ZLIB_AVAILABLE__)

include(ExternalProject)
ExternalProject_Add(HTSLIB
        PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
        STAMP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/build/stamp
        GIT_REPOSITORY https://github.com/samtools/htslib.git
        GIT_TAG "master"

        UPDATE_COMMAND ""

        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/samtools/htslib
        CONFIGURE_COMMAND autoheader && autoconf && ${CMAKE_CURRENT_SOURCE_DIR}/samtools/htslib/configure --prefix=${CMAKE_CURRENT_SOURCE_DIR}/samtools/htslib --disable-lzma
        BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/samtools/htslib
        BUILD_COMMAND make && make test
)

add_subdirectory(statgen)
set(LIBSTATGEN_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/statgen")
add_subdirectory(samtools)
add_subdirectory(libVcf)
#add_subdirectory(libEigen)


find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(SOURCE_FILES main.cpp SVDcalculator.cpp ContaminationEstimator.cpp MathGenMin.cpp MathGold.cpp Random.cpp SimplePileupViewer.cpp params.cpp )#SVDcalculator.cpp SVDcalculator.h)
include_directories(statgen ${CMAKE_CURRENT_SOURCE_DIR}/samtools/htslib/include/ samtools libVcf Eigen)
add_executable(VerifyBamID ${SOURCE_FILES})

find_library(ZLIB z HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT ZLIB)
    message(FATAL_ERROR "libz library not found")
endif()
find_library(CURLLIB curl HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT CURLLIB)
    message(FATAL_ERROR "libcurl library not found")
endif()
find_library(BZ2LIB bz2 HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT BZ2LIB)
    message(FATAL_ERROR "libbz2 library not found")
endif()
find_library(SSLLIB ssl HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT SSLLIB)
    message(FATAL_ERROR "libssl library not found")
endif()
find_library(CRYPTOLIB crypto HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT CRYPTOLIB)
    message(FATAL_ERROR "libcrypto library not found")
endif()

#add_dependencies(VerifyBamID HTSLIB)
target_link_libraries(VerifyBamID statgen Vcf samtools ${CMAKE_CURRENT_SOURCE_DIR}/samtools/htslib/libhts.a ${ZLIB} ${BZ2LIB} ${CURLLIB} ${SSLLIB} ${CRYPTOLIB})

enable_testing()
add_test( NAME myTest1
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
          COMMAND VerifyBamID --UDPath ./resource/test/hapmap_3.3.b37.dat.UD --BamFile ./resource/test/test.bam --BedPath resource/test/hapmap_3.3.b37.dat.bed --MeanPath ./resource/test/hapmap_3.3.b37.dat.mu --Reference ./resource/test/chr20.fa.gz )
add_test( NAME myTest2
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
          COMMAND sh ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/run.plot.sh -i ./resource/test/hapmap_3.3.b37.dat.V -o ./resource/test/hapmap -r 1000g -g grey)
