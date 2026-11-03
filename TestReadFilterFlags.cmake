# Exercise the actual command-line option binding and BAM pileup, including
# overlap correction. Keep every output under the build tree.
file(MAKE_DIRECTORY "${TEST_DIR}")
set(prefix "${TEST_DIR}/reads")
execute_process(COMMAND "${MAKE_BAM}" "${prefix}" RESULT_VARIABLE status)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "Could not generate read-filter fixture")
endif()
file(WRITE "${prefix}.bed" "chr1\t99\t100\tA\tC\n")
file(WRITE "${prefix}.mu" "chr1:100\t1\n")
file(WRITE "${prefix}.UD" "0\t0\n")

foreach(mode default zero proper_pair)
    set(options)
    set(expected_depth 2)
    if(mode STREQUAL "zero")
        set(options --incl-flags 0)
    elseif(mode STREQUAL "proper_pair")
        set(options --incl-flags 2)
        set(expected_depth 1)
    endif()
    execute_process(
        COMMAND "${VERIFY_BAM_ID}" --BamFile "${prefix}.bam"
            --Reference "${prefix}.fa" --SVDPrefix "${prefix}"
            --Output "${TEST_DIR}/${mode}" --OutputPileup
            --WithinAncestry --FixPC 0:0 --DisableSanityCheck ${options}
        RESULT_VARIABLE status OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
    if(NOT status EQUAL 0)
        message(FATAL_ERROR "${mode} failed: ${stdout}\n${stderr}")
    endif()
    file(STRINGS "${TEST_DIR}/${mode}.Pileup" rows)
    list(LENGTH rows row_count)
    if(NOT row_count EQUAL 1)
        message(FATAL_ERROR "${mode}: expected one marker, got ${rows}")
    endif()
    string(REPLACE "\t" ";" fields "${rows}")
    list(GET fields 3 depth)
    if(NOT depth EQUAL expected_depth)
        message(FATAL_ERROR "${mode}: expected depth ${expected_depth}, got ${depth}")
    endif()
endforeach()
