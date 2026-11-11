cmake_minimum_required(VERSION 3.10)
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

# statgen reads gzip files through zlib itself, whichever htslib is used.
find_package(ZLIB)
if (NOT ZLIB_FOUND)
    message(FATAL_ERROR "zlib development files not found")
endif ()

# HTS_INCLUDE_DIRS and HTS_LIBRARIES name an installed htslib. A -D on the
# command line leaves its cache entry untyped until find_path() or
# find_library() below adopts it, which tells a value given in this configure
# apart from one cached by an earlier system-htslib configure.
set(HTS_VARIABLES_GIVEN OFF)
foreach (variable HTS_INCLUDE_DIRS HTS_LIBRARIES)
    if (DEFINED ${variable})
        get_property(type CACHE ${variable} PROPERTY TYPE)
        if (type STREQUAL "UNINITIALIZED")
            set(HTS_VARIABLES_GIVEN ON)
        endif ()
    endif ()
endforeach ()
if (NOT DEFINED USE_SYSTEM_HTSLIB AND HTS_VARIABLES_GIVEN)
    set(USE_SYSTEM_HTSLIB_DEFAULT ON)
else ()
    set(USE_SYSTEM_HTSLIB_DEFAULT OFF)
endif ()
option(USE_SYSTEM_HTSLIB "Link an installed htslib instead of downloading and building htslib and libdeflate" ${USE_SYSTEM_HTSLIB_DEFAULT})

if (USE_SYSTEM_HTSLIB)
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

    set(HTSLIB_INCLUDE_DIRS ${HTS_INCLUDE_DIRS})
    set(HTSLIB_LINK_LIBRARIES ${HTS_LIBRARIES} ${ZLIB} ${BZIP2} ${LZMA} ${CURLLIB})
    message(STATUS "Using the installed htslib: ${HTS_LIBRARIES}")
else ()
    if (HTS_VARIABLES_GIVEN)
        message(FATAL_ERROR "HTS_INCLUDE_DIRS and HTS_LIBRARIES only apply with -DUSE_SYSTEM_HTSLIB=ON")
    endif ()
    if (HTS_INCLUDE_DIRS OR HTS_LIBRARIES)
        message(WARNING "Building the bundled htslib instead of the installed one cached by an "
            "earlier configure. Pass -DUSE_SYSTEM_HTSLIB=ON to keep using the installed htslib.")
        unset(HTS_INCLUDE_DIRS CACHE)
        unset(HTS_LIBRARIES CACHE)
    endif ()
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/BundledHtslib.cmake)
    message(STATUS "Building the bundled htslib ${HTSLIB_VERSION} and libdeflate ${LIBDEFLATE_VERSION}")
endif ()

# Every target that uses htslib goes through here, so each one compiles and
# links against the same htslib.
function(target_link_htslib target)
    target_include_directories(${target} PRIVATE ${HTSLIB_INCLUDE_DIRS})
    target_link_libraries(${target} ${HTSLIB_LINK_LIBRARIES})
    if (HTSLIB_BUILD_TARGET)
        add_dependencies(${target} ${HTSLIB_BUILD_TARGET})
    endif ()
endfunction()

add_subdirectory(statgen)
set(LIBSTATGEN_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/statgen")
add_subdirectory(samtools)
add_subdirectory(libVcf)
set(SOURCE_FILES main.cpp SVDcalculator.cpp ContaminationEstimator.cpp MathGenMin.cpp MathGold.cpp Random.cpp SimplePileupViewer.cpp MethylationModel.cpp ReadType.cpp params.cpp )
include_directories(statgen samtools libVcf Eigen)
add_executable(VerifyBamID ${SOURCE_FILES})

target_link_libraries(VerifyBamID statgen Vcf samtools)
target_link_htslib(VerifyBamID)

add_executable(TestReadVcf TestReadVcf.cpp SVDcalculator.cpp)
target_link_libraries(TestReadVcf statgen Vcf samtools)
target_link_htslib(TestReadVcf)

add_executable(TestGramSVD TestGramSVD.cpp SVDcalculator.cpp)
target_include_directories(TestGramSVD PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(TestGramSVD statgen Vcf samtools)
target_link_htslib(TestGramSVD)

add_executable(TestMethylationModel TestMethylationModel.cpp MethylationModel.cpp)
target_include_directories(TestMethylationModel PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_htslib(TestMethylationModel)

add_executable(TestReadType TestReadType.cpp ReadType.cpp)
target_link_libraries(TestReadType statgen)
target_link_htslib(TestReadType)

add_executable(TestHtslibFeatures TestHtslibFeatures.cpp)
target_link_htslib(TestHtslibFeatures)

# Test-only helper: synthesizes a directional methyl-seq BAM from an aligned BAM
# so the end-to-end --Methylation test can generate its input programmatically.
add_executable(MethylateTestBam MethylateTestBam.cpp)
target_link_htslib(MethylateTestBam)

add_executable(TestMixLikelihood TestMixLikelihood.cpp ContaminationEstimator.cpp
        MathGenMin.cpp MathGold.cpp Random.cpp SimplePileupViewer.cpp MethylationModel.cpp)
target_link_libraries(TestMixLikelihood statgen Vcf samtools)
target_link_htslib(TestMixLikelihood)

enable_testing()
add_executable(MakeMaxDepthTestBam MakeMaxDepthTestBam.cpp)
target_link_htslib(MakeMaxDepthTestBam)
foreach(mode caps_every_marker zero_means_no_limit negative_is_rejected)
    add_test(NAME testMaxDepth_${mode}
            COMMAND ${CMAKE_COMMAND} -DMODE=${mode}
                -DVERIFY_BAM_ID=$<TARGET_FILE:VerifyBamID>
                -DMAKE_BAM=$<TARGET_FILE:MakeMaxDepthTestBam>
                -DTEST_DIR=${CMAKE_CURRENT_BINARY_DIR}/max-depth-test-${mode}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/TestMaxDepth.cmake)
endforeach()
add_executable(MakeReadFilterTestBam MakeReadFilterTestBam.cpp)
target_link_htslib(MakeReadFilterTestBam)
add_test(NAME testReadFilterFlags
        COMMAND ${CMAKE_COMMAND}
            -DVERIFY_BAM_ID=$<TARGET_FILE:VerifyBamID>
            -DMAKE_BAM=$<TARGET_FILE:MakeReadFilterTestBam>
            -DTEST_DIR=${CMAKE_CURRENT_BINARY_DIR}/read-filter-test
            -P ${CMAKE_CURRENT_SOURCE_DIR}/TestReadFilterFlags.cmake)

# Long-read pileup tests share data generated by the setup test.
add_executable(MakeLongReadTestData MakeLongReadTestData.cpp)
target_link_htslib(MakeLongReadTestData)
set(LONG_READ_TEST_ARGS
        -DVERIFY_BAM_ID=$<TARGET_FILE:VerifyBamID>
        -DMAKE_DATA=$<TARGET_FILE:MakeLongReadTestData>
        -DTEST_DIR=${CMAKE_CURRENT_BINARY_DIR}/long-read-test
        -P ${CMAKE_CURRENT_SOURCE_DIR}/TestLongReadPileup.cmake)
add_test(NAME testLongReadPileupSetup COMMAND ${CMAKE_COMMAND} -DMODE=setup ${LONG_READ_TEST_ARGS})
add_test(NAME testLongReadPileupFromBam COMMAND ${CMAKE_COMMAND} -DMODE=bam ${LONG_READ_TEST_ARGS})
add_test(NAME testLongReadPileupFromCram COMMAND ${CMAKE_COMMAND} -DMODE=cram ${LONG_READ_TEST_ARGS})
add_test(NAME testReadTypeShortKeepsBaqAndAdjustMq COMMAND ${CMAKE_COMMAND} -DMODE=short_read_type ${LONG_READ_TEST_ARGS})
add_test(NAME testLongReadsKeepExplicitAdjustMq COMMAND ${CMAKE_COMMAND} -DMODE=explicit_adjust_mq ${LONG_READ_TEST_ARGS})
add_test(NAME testUnknownReadTypeIsRejected COMMAND ${CMAKE_COMMAND} -DMODE=unknown_read_type ${LONG_READ_TEST_ARGS})
add_test(NAME testAdjustMqBelowMinusOneIsRejected COMMAND ${CMAKE_COMMAND} -DMODE=adjust_mq_below_minus_one ${LONG_READ_TEST_ARGS})
add_test(NAME testTruncatedBamStopsTheRun COMMAND ${CMAKE_COMMAND} -DMODE=truncated_bam ${LONG_READ_TEST_ARGS})
add_test(NAME testLongReadsKeepBasesNearIndelsWhenDiscardIsOff
        COMMAND ${CMAKE_COMMAND} -DMODE=discard_near_indel_off ${LONG_READ_TEST_ARGS})
add_test(NAME testShortReadsCountEveryAlignment
        COMMAND ${CMAKE_COMMAND} -DMODE=short_reads_count_every_alignment ${LONG_READ_TEST_ARGS})
add_test(NAME testShortReadsKeepExplicitDiscardNearIndel
        COMMAND ${CMAKE_COMMAND} -DMODE=explicit_discard_near_indel ${LONG_READ_TEST_ARGS})
add_test(NAME testDiscardNearIndelBelowMinusOneIsRejected
        COMMAND ${CMAKE_COMMAND} -DMODE=discard_near_indel_below_minus_one ${LONG_READ_TEST_ARGS})
add_test(NAME testLargestDiscardNearIndelSparesReadsWithoutIndels
        COMMAND ${CMAKE_COMMAND} -DMODE=huge_discard_near_indel ${LONG_READ_TEST_ARGS})
add_test(NAME testDiscardNearIndelWarnsWithPileupFile
        COMMAND ${CMAKE_COMMAND} -DMODE=discard_near_indel_with_pileup_file ${LONG_READ_TEST_ARGS})
set_tests_properties(testLongReadPileupFromBam testLongReadPileupFromCram
        testReadTypeShortKeepsBaqAndAdjustMq testLongReadsKeepExplicitAdjustMq testUnknownReadTypeIsRejected
        testAdjustMqBelowMinusOneIsRejected testTruncatedBamStopsTheRun
        testLongReadsKeepBasesNearIndelsWhenDiscardIsOff testShortReadsCountEveryAlignment
        testShortReadsKeepExplicitDiscardNearIndel testDiscardNearIndelBelowMinusOneIsRejected
        testLargestDiscardNearIndelSparesReadsWithoutIndels testDiscardNearIndelWarnsWithPileupFile
        PROPERTIES DEPENDS testLongReadPileupSetup)

add_test(NAME testMixLikelihood COMMAND TestMixLikelihood)
add_test(NAME testReadVcf
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND TestReadVcf resource/test/test_readvcf.vcf)
add_test(NAME testGramSVD
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND TestGramSVD)
add_test(NAME testMethylationModel
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND TestMethylationModel)
add_test(NAME testReadType
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND TestReadType)
# Whichever htslib is linked must read CRAM files compressed with lzma or
# bzip2; the bundled one must also open remote files and use libdeflate.
add_test(NAME testHtslibSupportsLzmaAndBzip2 COMMAND TestHtslibFeatures lzma bzip2)
if (NOT USE_SYSTEM_HTSLIB)
    add_test(NAME testBundledHtslibHasRemoteAccessAndLibdeflate
            COMMAND TestHtslibFeatures libcurl s3 gcs libdeflate)
endif ()
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

# End-to-end --Methylation test. The methyl-seq input is synthesized from the
# existing test.bam at test time (no committed test data): MethylateTestBam
# applies the directional bisulfite conversion, then VerifyBamID --Methylation
# is diffed against a golden Ancestry file. Chained via DEPENDS.
add_test(NAME methylMakeBam
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND MethylateTestBam resource/test/test.bam resource/test/test.methyl.bam)
add_test(NAME testMethylationRun
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND VerifyBamID --DisableSanityCheck --Methylation --BamFile resource/test/test.methyl.bam --SVDPrefix resource/test/hapmap_3.3.b37.dat --Reference resource/test/chr20.fa.gz --NumPC 2 --Output result.Methylation)
set_tests_properties(testMethylationRun PROPERTIES DEPENDS methylMakeBam)
add_test(NAME testMethylationDiff
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND sh -c "diff resource/test/expected/result.methyl.Ancestry result.Methylation.Ancestry && rm result.Methylation.Ancestry result.Methylation.selfSM resource/test/test.methyl.bam resource/test/test.methyl.bam.bai")
set_tests_properties(testMethylationDiff PROPERTIES DEPENDS testMethylationRun)

# --ignoreRG regression (issue #50). resource/test/test.twoRG.bam carries two
# read groups with distinct SM tags: by default the tool rejects it as
# multi-sample, while --ignoreRG wires MPLP_IGNORE_RG so all reads collapse
# into one sample group named by the input file.
add_test(NAME testIgnoreRGDefaultFails
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND VerifyBamID --DisableSanityCheck --BamFile resource/test/test.twoRG.bam --SVDPrefix resource/test/hapmap_3.3.b37.dat --Reference resource/test/chr20.fa.gz --NumPC 2 --Output result.IgnoreRGDefault)
set_tests_properties(testIgnoreRGDefaultFails PROPERTIES WILL_FAIL TRUE)
add_test(NAME testIgnoreRGRun
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND VerifyBamID --DisableSanityCheck --ignoreRG --BamFile resource/test/test.twoRG.bam --SVDPrefix resource/test/hapmap_3.3.b37.dat --Reference resource/test/chr20.fa.gz --NumPC 2 --Output result.IgnoreRG)
add_test(NAME testIgnoreRGDiff
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMAND sh -c "test `tail -n +2 result.IgnoreRG.selfSM | wc -l` -eq 1 && test \"x`head -2 result.IgnoreRG.selfSM | tail -1 | cut -f1`\" = \"xresource/test/test.twoRG.bam\" && rm -f result.IgnoreRG.Ancestry result.IgnoreRG.selfSM")
set_tests_properties(testIgnoreRGDiff PROPERTIES DEPENDS testIgnoreRGRun)

#add_test( NAME myTestPlot
#          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#          COMMAND sh ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/run.plot.sh -i ${CMAKE_CURRENT_SOURCE_DIR}/resource/test/hapmap_3.3.b37.dat.V -o ${CMAKE_CURRENT_SOURCE_DIR}/resource/test/hapmap -r 1000g -g grey)
