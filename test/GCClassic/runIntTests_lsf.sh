#!/bin/bash

#### Liam: add LSF tags here!

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: runIntTests_slurm.sh
#
# !DESCRIPTION: Runs integration tests on the various GEOS-Chem Classic
#  run directories (using the SLURM scheduler).
#\\
#\\
# !CALLING SEQUENCE:
#  bsub runIntTests_slurm.sh
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

# Load common functions
. ${ROOT}/commonFunctionsForTests.sh

# Count the number of tests to be done = number of run directories
NUM_TESTS=$(count_rundirs ${ROOT})

#============================================================================
# Initialize results logfile
#============================================================================

# Results logfile name
RESULTS=${ROOT}/logs/results.log

# Print header to results log file
print_to_log "${LINE}"                                       ${RESULTS}
print_to_log "GEOS-Chem Integration Test Results"            ${RESULTS}
print_to_log ""                                              ${RESULTS}
print_to_log "Using ${OMP_NUM_THREADS} OpenMP threads"       ${RESULTS}
print_to_log "Number of tests to be performed: ${NUM_TESTS}" ${RESULTS}
print_to_log "${LINE}"                                       ${RESULTS}

#============================================================================
# Configure and compile code in each GEOS_Chem run directory
#============================================================================
print_to_log " "                   ${RESULTS}
print_to_log "Compiliation tests:" ${RESULTS}
print_to_log "${LINELC}"           ${RESULTS}

# Loop over rundirs and compile code
for RUNDIR in *; do
    if [[ -d ${RUNDIR} && "x${RUNDIR}" != "xlogs" ]]; then
	LOG=${ROOT}/logs/compile.${RUNDIR}.log
	config_and_build ${ROOT} ${RUNDIR} ${LOG} ${RESULTS}
    fi
done

#============================================================================
# Run the GEOS-Chem executable in each GEOS-Chem run directory
#============================================================================
print_to_log " "                 ${RESULTS}
print_to_log "Execution tests:"  ${RESULTS}
print_to_log "${LINELC}"         ${RESULTS}

# Loop over rundirs and run GEOS-Chem
for RUNDIR in *; do
    if [[ -d ${RUNDIR} && "x${RUNDIR}" != "xlogs" ]]; then
	LOG=${ROOT}/logs/run.${RUNDIR}.log
	run_gcclassic ${ROOT} ${RUNDIR} ${LOG} ${RESULTS}
    fi
done

#============================================================================
# Check the number of simulations that have passed
#============================================================================

# Look for matches
CMP_PASSED=$(count_matches_in_log "${CMP_PASS_STR}" ${RESULTS})
CMP_FAILED=$(count_matches_in_log "${CMP_FAIL_STR}" ${RESULTS})
RUN_PASSED=$(count_matches_in_log "${RUN_PASS_STR}" ${RESULTS})
RUN_FAILED=$(count_matches_in_log "${RUN_FAIL_STR}" ${RESULTS})

# Print summary to log
print_to_log " "                                        ${RESULTS}
print_to_log "Summary of test results:"                 ${RESULTS}
print_to_log "${LINELC}"                                ${RESULTS}
print_to_log "Complilation tests passed: ${CMP_PASSED}" ${RESULTS}
print_to_log "Complilation tests failed: ${CMP_FAILED}" ${RESULTS}
print_to_log "Execution    tests passed: ${RUN_PASSED}" ${RESULTS}
print_to_log "Execution    tests failed: ${RUN_FAILED}" ${RESULTS}

# Check if all tests passed
if [[ "x${CMP_PASSED}" == "x${NUM_TESTS}" ]] && \
   [[ "x${RUN_PASSED}" == "x${NUM_TESTS}" ]]; then
    print_to_log ""                             ${RESULTS}
    print_to_log "%%%%%%%%%%%%%%%%%%%%%%%%%%%%" ${RESULTS}
    print_to_log "%%%  All tests passed!   %%%" ${RESULTS}
    print_to_log "%%%%%%%%%%%%%%%%%%%%%%%%%%%%" ${RESULTS}
fi

#============================================================================
# Cleanup and quit
#============================================================================

# Free local variables
unset LOG
unset NUM_TESTS
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
