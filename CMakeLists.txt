cmake_minimum_required(VERSION 2.8.12)
project(VerifyBamID)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3 -fPIC -Wno-format-security")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -O3  -fPIC -Wno-unused-variable -Wno-unused-result ")

add_definitions(-D__STDC_LIMIT_MACROS)
add_definitions(-D__ZLIB_AVAILABLE__)

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

find_library(CURLLIB curl HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ /usr/local/lib/)
if (NOT CURLLIB)
    message(FATAL_ERROR "libcurl library not found")
endif ()

find_library(ZLIB z HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ /usr/local/lib/ /usr/local/Cellar/zlib/1.2.11/lib/)
if (NOT ZLIB)
    message(FATAL_ERROR "libz library not found")
endif ()

find_path(HTS_INCLUDE_DIRS htslib/bgzf.h HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ /usr/local/include/)
if (NOT HTS_INCLUDE_DIRS)
    message(FATAL_ERROR "libhts HTS_INCLUDE_DIRS not found")
endif ()

find_library(HTS_LIBRARIES hts HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ /usr/local/lib/)
if (NOT HTS_LIBRARIES)
    message(FATAL_ERROR "libhts HTS_LIBRARIES not found")
endif ()

find_package(BZip2)
if (NOT BZIP2_FOUND)
    message(FATAL_ERROR "libbz2 library not found")
else ()
    include_directories(${BZIP2_INCLUDE_DIRS})
    set(BZIP2 ${BZIP2_LIBRARIES})
endif ()

find_library(LZMA lzma HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ /usr/local/lib/)
if(NOT LZMA)
    message(WARNING "liblzma library not found, if you specified --disable-lzma when compiling libhts, please ignore this warning")
    set(LZMA "")
else()
    set(LZMA ${LZMA_LIBRARIES})
endif()


add_subdirectory(statgen)
set(LIBSTATGEN_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/statgen")
add_subdirectory(samtools)
add_subdirectory(libVcf)
set(SOURCE_FILES main.cpp SVDcalculator.cpp ContaminationEstimator.cpp MathGenMin.cpp MathGold.cpp Random.cpp SimplePileupViewer.cpp params.cpp )
include_directories(statgen ${HTS_INCLUDE_DIRS} samtools libVcf Eigen)
add_executable(VerifyBamID ${SOURCE_FILES})

target_link_libraries(VerifyBamID statgen Vcf ${HTS_LIBRARIES} samtools ${ZLIB} ${BZIP2} ${LZMA} ${CURLLIB})

enable_testing()
add_test(NAME myTest1
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND VerifyBamID --DisableSanityCheck --BamFile resource/test/test.bam --SVDPrefix resource/test/hapmap_3.3.b37.dat --Reference resource/test/chr20.fa.gz --NumPC 2)
#add_test( NAME myTest2
#          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#          COMMAND sh ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/run.plot.sh -i ${CMAKE_CURRENT_SOURCE_DIR}/resource/test/hapmap_3.3.b37.dat.V -o ${CMAKE_CURRENT_SOURCE_DIR}/resource/test/hapmap -r 1000g -g grey)
