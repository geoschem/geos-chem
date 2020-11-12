#!/bin/bash

#SBATCH -c 24
#SBATCH -N 1
#SBATCH -t 0-18:00
#SBATCH -p huce_cascade
#SBATCH --mem=60000
#SBATCH --mail-type=END

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: intTestExecute_slurm.sh
#
# !DESCRIPTION: Runs execution tests on various GEOS-Chem Classic
#  run directories (using the SLURM scheduler).
#\\
#\\
# !CALLING SEQUENCE:
#  sbatch intTestExecute_slurm.sh
#
# !REVISION HISTORY:
#  03 Nov 2020 - R. Yantosca - Initial version
#  See the subsequent Git history with the gitk browser!
#EOP
#------------------------------------------------------------------------------
#BOC

#============================================================================
# Global variable and function definitions
#============================================================================

# Get the long path of this folder
ROOT=`pwd -P`

# Load base computational environment
. ~/.bashrc

# Load software modules and OpenMP settings
. ${ROOT}/gcclassic_env.sh

# Load common functions for tests
. ${ROOT}/commonFunctionsForTests.sh

# Count the number of tests to be done = number of run directories
NUM_TESTS=$(count_rundirs ${ROOT})

#============================================================================
# Initialize results logfile
#============================================================================

# Results logfile name
RESULTS="${ROOT}/logs/results.execute.log"
rm -f ${RESULTS}

# Print header to results log file
print_to_log "${SEP_MAJOR}"                              ${RESULTS}
print_to_log "GEOS-Chem Classic: Execution Test Results" ${RESULTS}
print_to_log ""                                          ${RESULTS}
print_to_log "Using ${OMP_NUM_THREADS} OpenMP threads"   ${RESULTS}
print_to_log "Number of execution tests: ${NUM_TESTS}"   ${RESULTS}
print_to_log "${SEP_MAJOR}"                              ${RESULTS}

#============================================================================
# Run the GEOS-Chem executable in each GEOS-Chem run directory
#============================================================================
print_to_log " "                 ${RESULTS}
print_to_log "Execution tests:"  ${RESULTS}
print_to_log "${SEP_MINOR}"      ${RESULTS}

# Loop over rundirs and run GEOS-Chem
# Keep track of the number of tests that passed & failed
let PASSED=0
let FAILED=0
let REMAIN=${NUM_TESTS}
for RUNDIR in *; do

    # Do the following if for only valid GEOS-Chem run dirs
    EXPR=$(is_valid_rundir "${ROOT}/${RUNDIR}")
    if [[ "x${EXPR}" == "xTRUE" ]]; then

	# Define log file
	LOG="${ROOT}/logs/execute.${RUNDIR}.log"
	rm -f ${LOG}

	# Messages for execution pass & fail
	PASSMSG="$RUNDIR${FILL:${#RUNDIR}}.....${EXE_PASS_STR}"
	FAILMSG="$RUNDIR${FILL:${#RUNDIR}}.....${EXE_FAIL_STR}"

	# Test if executable is present
	if [[ -f ${ROOT}/${RUNDIR}/gcclassic ]]; then

	    # Change to this run directory; remove leftover log file
	    cd ${ROOT}/${RUNDIR}

	    # Run the code if the executable is present.  Then update the
	    # pass/fail counters and write a message to the results log file.
	    ./gcclassic >> ${LOG} 2>&1
	    if [[ $? -eq 0 ]]; then
		let PASSED++
		if [[ "x${RESULTS}" != "x" ]]; then
		    print_to_log "${PASSMSG}" ${RESULTS}
		fi
	    else
		let FAILED++
		if [[ "x${RESULTS}" != "x" ]]; then
		    print_to_log "${FAILMSG}" ${RESULTS}
		fi
	    fi

	    # Change to root directory for next iteration
	    cd ${ROOT}

	else

	    # If the executable is missing, update the "fail" counter
	    # and write the "failed" message to the results log file.
	    let FAILED++
	    if [[ "x${RESULTS}" != "x" ]]; then
		print_to_log "${FAILMSG}" ${results}
	    fi
	fi

	# Decrement the count of remaining tests
	let REMAIN--
    fi
done

#============================================================================
# Check the number of simulations that have passed
#============================================================================

# Print summary to log
print_to_log " "                                            ${RESULTS}
print_to_log "Summary of test results:"                     ${RESULTS}
print_to_log "${SEP_MINOR}"                                 ${RESULTS}
print_to_log "Execution tests passed: ${PASSED}"            ${RESULTS}
print_to_log "Execution tests failed: ${FAILED}"            ${RESULTS}
print_to_log "Execution tests not yet completed: ${REMAIN}" ${RESULTS}

# Check if all tests passed
if [[ "x${PASSED}" == "x${NUM_TESTS}" ]]; then
    print_to_log ""                                         ${RESULTS}
    print_to_log "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"    ${RESULTS}
    print_to_log "%%%  All execution tests passed!  %%%"    ${RESULTS}
    print_to_log "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"    ${RESULTS}
fi

#============================================================================
# Cleanup and quit
#============================================================================

# Free local variables
unset FAILED
unset FAILMSG
unset LOG
unset NUM_TESTS
unset PASSED
unset PASSMSG
unset REMAIN
unset RESULTS
unset ROOT

# Free imported global variables
unset FILL
unset LINE
unset LINELC
unset SED_INPUT_GEOS_1
unset SED_INPUT_GEOS_2
unset SED_HISTORY_RC
unset CMP_PASS_STR
unset CMP_FAIL_STR
unset RUN_PASS_STR
unset RUN_FAIL_STR
#EOC
