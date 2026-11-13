# Checks --max-depth on a BAM from MakeMaxDepthTestBam, which has five reads at
# each of two markers. MODE is one of:
#   caps_every_marker     --max-depth 3 keeps three of the five reads at both markers
#   zero_means_no_limit   --max-depth 0 keeps all five reads at both markers
#   negative_is_rejected  --max-depth -1 is an error
# Keep every output under the build tree.
file(REMOVE_RECURSE "${TEST_DIR}")
file(MAKE_DIRECTORY "${TEST_DIR}")
set(prefix "${TEST_DIR}/reads")
execute_process(COMMAND "${MAKE_BAM}" "${prefix}" RESULT_VARIABLE status)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "Could not generate max-depth fixture")
endif()
file(WRITE "${prefix}.bed" "chr1\t499\t500\tA\tC\nchr1\t1499\t1500\tA\tC\n")
file(WRITE "${prefix}.mu" "chr1:500\t1\nchr1:1500\t1\n")
file(WRITE "${prefix}.UD" "0\t0\n0\t0\n")

# Runs VerifyBamID with `max_depth` and sets `status` and `log` in the caller.
function(run_verify_bam_id max_depth)
    execute_process(
        COMMAND "${VERIFY_BAM_ID}" --BamFile "${prefix}.bam"
            --Reference "${prefix}.fa" --SVDPrefix "${prefix}"
            --Output "${TEST_DIR}/out" --OutputPileup --max-depth ${max_depth}
            --WithinAncestry --FixPC 0:0 --DisableSanityCheck
        RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
    set(status "${result}" PARENT_SCOPE)
    set(log "${stdout}${stderr}" PARENT_SCOPE)
endfunction()

# Sets `depths` in the caller to the depth column of each .Pileup row.
function(read_depths)
    file(STRINGS "${TEST_DIR}/out.Pileup" rows)
    set(result)
    foreach(row ${rows})
        string(REPLACE "\t" ";" fields "${row}")
        list(GET fields 3 depth)
        list(APPEND result ${depth})
    endforeach()
    set(depths "${result}" PARENT_SCOPE)
endfunction()

if(MODE STREQUAL "caps_every_marker")
    run_verify_bam_id(3)
    if(NOT status EQUAL 0)
        message(FATAL_ERROR "--max-depth 3 failed:\n${log}")
    endif()
    read_depths()
    if(NOT depths STREQUAL "3;3")
        message(FATAL_ERROR "--max-depth 3: expected three reads at both markers, got ${depths}")
    endif()
elseif(MODE STREQUAL "zero_means_no_limit")
    run_verify_bam_id(0)
    if(NOT status EQUAL 0)
        message(FATAL_ERROR "--max-depth 0 failed:\n${log}")
    endif()
    read_depths()
    if(NOT depths STREQUAL "5;5")
        message(FATAL_ERROR "--max-depth 0: expected all five reads at both markers, got ${depths}")
    endif()
elseif(MODE STREQUAL "negative_is_rejected")
    run_verify_bam_id(-1)
    string(FIND "${log}" "--max-depth must be 0 (no limit) or positive" found)
    if(status EQUAL 0 OR found EQUAL -1)
        message(FATAL_ERROR "--max-depth -1: expected an error, got status ${status}:\n${log}")
    endif()
else()
    message(FATAL_ERROR "Unknown MODE '${MODE}'")
endif()
