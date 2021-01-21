#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: intTestResults.sh
#
# !DESCRIPTION: Displays results from the execution phase of the
#  GCHP integration tests.
#\\
#\\
# !CALLING SEQUENCE:
#  sbatch intTestResults.sh
#
# !REVISION HISTORY:
#  03 Nov 2020 - R. Yantosca - Initial version
#  See the subsequent Git history with the gitk browser!
#EOP
#------------------------------------------------------------------------------
#BOC

#============================================================================
# Variable and function definitions
#============================================================================

# Get the long path of this folder
root=$(pwd -P)

# Load common functions for tests
. ${root}/commonFunctionsForTests.sh

# Number of execution tests = number of run directories
numTests=$(count_rundirs ${root})

# Results logfile name
results="${root}/logs/results.execute.log"
rm -f ${results}

#============================================================================
# Initialize results logfile
#============================================================================

# Print header to results log file
print_to_log "${SEP_MAJOR}"                             ${results}
print_to_log "GCHP: Execution Test Results"          ${results}
print_to_log ""                                         ${results}
print_to_log "Number of execution tests: ${numTests}"   ${results}
print_to_log "${SEP_MAJOR}"                             ${results}

#============================================================================
# Configure and compile code in each GEOS_Chem run directory
#============================================================================
print_to_log " "                ${results}
print_to_log "Execution tests:" ${results}
print_to_log "${SEP_MINOR}"     ${results}

# Change to the top-level build directory
cd ${root}

# Keep track of the number of tests that passed & failed
let passed=0
let failed=0
let remain=${numTests}

# Get a list of GCHP run directories
dirlist=$(ls -1 | grep gchp_)

# Loop over all of the execution logs
for dir in ${dirlist}; do

    # Log file path
    path="execute.${dir}.log"

    # If the log file path exists, then get the filename part
    file=$(basename ${path})
    runDir="${file%.*}"
    runDir="${runDir#*.}"

    # Create sucess and failure messages
    passMsg="$runDir${FILL:${#runDir}}.....${EXE_PASS_STR}"
    failMsg="$runDir${FILL:${#runDir}}.....${EXE_FAIL_STR}"
    naMsg="$runDir${FILL:${#runDir}}.....Not Completed"

    # Test if the log file exists
    if [[ -f ${path} ]]; then

	# Look for the text ----EXTDATA, which shows up
	# at the end of a successful GCHP job
	grep -e ----EXTDATA ${path} > /dev/null
	if [[ $? -eq 0 ]]; then
	    let passed++
	    print_to_log "${passMsg}" ${results}
	else
    	    let failed++
	    print_to_log "${failMsg}" ${results}
	fi
	let remain--

    else

	# If the log file path does not exist,
	# then indicate the test is still to be done
	print_to_log "${naMsg}" ${results}

    fi
done

#============================================================================
# Check the number of simulations that have passed
#============================================================================

# Print summary to log
print_to_log " "                                          ${results}
print_to_log "Summary of execution test results:"         ${results}
print_to_log "${SEP_MINOR}"                               ${results}
print_to_log "Execution tests passed:        ${passed}"   ${results}
print_to_log "Execution tests failed:        ${failed}"   ${results}
print_to_log "Execution tests not completed: ${remain}"   ${results}

# Check if all tests passed
if [[ "x${passed}" == "x${numTests}" ]]; then
    print_to_log ""                                      ${results}
    print_to_log "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" ${results}
    print_to_log "%%%  All execution tests passed!  %%%" ${results}
    print_to_log "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" ${results}
fi

#============================================================================
# Cleanup and quit
#============================================================================

# Free local variables
unset failed
unset file
unset runDir
unset log
unset numTests
unset options
unset passed
unset path
unset remain
unset results
unset root

# Free imported variables
unset FILL
unset SEP_MAJOR
unset SEP_MINOR
unset SED_INPUT_GEOS_1
unset SED_INPUT_GEOS_2
unset SED_HISTORY_RC
unset CMP_PASS_STR
unset CMP_FAIL_STR
unset EXE_PASS_STR
unset EXE_FAIL_STR
#EOC
