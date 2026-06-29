cmake_minimum_required(VERSION 3.5)
project(VerifyBamID)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3 -fPIC -Wno-format-security")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -O3 -fPIC -Wno-unused-variable -Wno-unused-result ")

# Set DEBUG compile definition conditionally based on build type
if(CMAKE_BUILD_TYPE MATCHES Debug)
    add_definitions(-DDEBUG=1)
else()
    add_definitions(-DDEBUG=0)
endif()

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

add_executable(TestReadVcf TestReadVcf.cpp SVDcalculator.cpp)
target_link_libraries(TestReadVcf statgen Vcf ${HTS_LIBRARIES} samtools ${ZLIB} ${BZIP2} ${LZMA} ${CURLLIB})

add_executable(TestGramSVD TestGramSVD.cpp SVDcalculator.cpp)
target_include_directories(TestGramSVD PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(TestGramSVD statgen Vcf ${HTS_LIBRARIES} samtools ${ZLIB} ${BZIP2} ${LZMA} ${CURLLIB})

enable_testing()
add_test(NAME testReadVcf
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND TestReadVcf resource/test/test_readvcf.vcf)
add_test(NAME testGramSVD
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND TestGramSVD)
add_test(NAME myTest1
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND VerifyBamID --DisableSanityCheck --BamFile resource/test/test.bam --SVDPrefix resource/test/hapmap_3.3.b37.dat --Reference resource/test/chr20.fa.gz --NumPC 2 --Output result.myTest1)
add_test(NAME myTest2
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND sh -c "diff resource/test/expected/result.Ancestry result.myTest1.Ancestry && rm result.myTest1.Ancestry result.myTest1.selfSM")
set_tests_properties(myTest2 PROPERTIES DEPENDS myTest1)
add_test(NAME myTest3
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND VerifyBamID --DisableSanityCheck --PileupFile resource/test/expected/result.Pileup --SVDPrefix resource/test/hapmap_3.3.b37.dat --Reference resource/test/chr20.fa.gz --NumPC 2 --Output result.myTest3)
add_test(NAME myTest4
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND sh -c "diff resource/test/expected/result.Ancestry result.myTest3.Ancestry && rm result.myTest3.Ancestry result.myTest3.selfSM")
set_tests_properties(myTest4 PROPERTIES DEPENDS myTest3)
add_test(NAME myTest5
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND VerifyBamID --DisableSanityCheck --PileupFile resource/test/test.LongRead.pileup --SVDPrefix resource/test/hapmap_3.3.b37.dat --Reference resource/test/chr20.fa.gz --NumPC 2 --Output result.myTest5)
add_test(NAME myTest6
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND sh -c "diff resource/test/expected/test.LongRead.pileup.Ancestry result.myTest5.Ancestry && rm result.myTest5.Ancestry result.myTest5.selfSM")
set_tests_properties(myTest6 PROPERTIES DEPENDS myTest5)

# Each new test pair below uses its own --Output prefix so tests can run in
# parallel without racing on the shared "result.*" filenames, and the Diff
# test is marked as DEPENDS on its Run so CTest always sequences them
# correctly regardless of scheduling.

# Test WithinAncestry model (OptimizeHomo)
add_test(NAME testWithinAncestryRun
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND VerifyBamID --DisableSanityCheck --BamFile resource/test/test.bam --SVDPrefix resource/test/hapmap_3.3.b37.dat --Reference resource/test/chr20.fa.gz --NumPC 2 --WithinAncestry --Output result.WithinAncestry)
add_test(NAME testWithinAncestryDiff
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND sh -c "diff resource/test/expected/test.WithinAncestry.Ancestry result.WithinAncestry.Ancestry && rm result.WithinAncestry.Ancestry result.WithinAncestry.selfSM")
set_tests_properties(testWithinAncestryDiff PROPERTIES DEPENDS testWithinAncestryRun)

# Test WithinAncestry + FixPC (OptimizeHomoFixedPC)
add_test(NAME testWithinAncestryFixPCRun
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND VerifyBamID --DisableSanityCheck --BamFile resource/test/test.bam --SVDPrefix resource/test/hapmap_3.3.b37.dat --Reference resource/test/chr20.fa.gz --NumPC 2 --WithinAncestry --FixPC 0.034756:0.0193 --Output result.WithinAncestryFixPC)
add_test(NAME testWithinAncestryFixPCDiff
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND sh -c "diff resource/test/expected/test.WithinAncestry.FixPC.Ancestry result.WithinAncestryFixPC.Ancestry && rm result.WithinAncestryFixPC.Ancestry result.WithinAncestryFixPC.selfSM")
set_tests_properties(testWithinAncestryFixPCDiff PROPERTIES DEPENDS testWithinAncestryFixPCRun)

# Test FixAlpha (OptimizeHomoFixedAlpha initial + OptimizeHeterFixedAlpha)
add_test(NAME testFixAlphaRun
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND VerifyBamID --DisableSanityCheck --BamFile resource/test/test.bam --SVDPrefix resource/test/hapmap_3.3.b37.dat --Reference resource/test/chr20.fa.gz --NumPC 2 --FixAlpha 0.1 --Output result.FixAlpha)
add_test(NAME testFixAlphaDiff
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND sh -c "diff resource/test/expected/test.FixAlpha.Ancestry result.FixAlpha.Ancestry && rm result.FixAlpha.Ancestry result.FixAlpha.selfSM")
set_tests_properties(testFixAlphaDiff PROPERTIES DEPENDS testFixAlphaRun)

# Test FixPC in BetweenAncestry mode (OptimizeHeterFixedPC)
add_test(NAME testHeterFixPCRun
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND VerifyBamID --DisableSanityCheck --BamFile resource/test/test.bam --SVDPrefix resource/test/hapmap_3.3.b37.dat --Reference resource/test/chr20.fa.gz --NumPC 2 --FixPC 0.034756:0.0193 --Output result.HeterFixPC)
add_test(NAME testHeterFixPCDiff
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND sh -c "diff resource/test/expected/test.HeterFixPC.Ancestry result.HeterFixPC.Ancestry && rm result.HeterFixPC.Ancestry result.HeterFixPC.selfSM")
set_tests_properties(testHeterFixPCDiff PROPERTIES DEPENDS testHeterFixPCRun)

#add_test( NAME myTestPlot
#          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#          COMMAND sh ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/run.plot.sh -i ${CMAKE_CURRENT_SOURCE_DIR}/resource/test/hapmap_3.3.b37.dat.V -o ${CMAKE_CURRENT_SOURCE_DIR}/resource/test/hapmap -r 1000g -g grey)
