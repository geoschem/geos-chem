#!/bin/bash

#SBATCH -c 24
#SBATCH -N 1
#SBATCH -t 0-12:00
#SBATCH -p huce_cascade
#SBATCH --mem=8000
#SBATCH --mail-type=END

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: intTestCompile_slurm.sh
#
# !DESCRIPTION: Runs compilation tests on various GEOS-Chem Classic
#  run directories (using the SLURM scheduler).
#\\
#\\
# !CALLING SEQUENCE:
#  sbatch intTestCompile_slurm.sh
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
RESULTS="${ROOT}/logs/results.compile.log"
rm -f ${RESULTS}

# Print header to results log file
print_to_log "${SEP_MAJOR}"                                ${RESULTS}
print_to_log "GEOS-Chem Classic: Compilation Test Results" ${RESULTS}
print_to_log ""                                            ${RESULTS}
print_to_log "Using ${OMP_NUM_THREADS} OpenMP threads"     ${RESULTS}
print_to_log "Number of compilation tests: ${NUM_TESTS}"   ${RESULTS}
print_to_log "${SEP_MAJOR}"                                ${RESULTS}

#============================================================================
# Configure and compile code in each GEOS_Chem run directory
#============================================================================
print_to_log " "                   ${RESULTS}
print_to_log "Compiliation tests:" ${RESULTS}
print_to_log "${SEP_MINOR}"        ${RESULTS}

# Loop over rundirs and compile code
# Keep track of the number of tests that passed & failed
let PASSED=0
let FAILED=0
let REMAIN=${NUM_TESTS}
for RUNDIR in *; do

    # Do the following only for valid GEOS-Chem run dirs
    EXPR=$(is_valid_rundir "${ROOT}/${RUNDIR}")
    if [[ "x${EXPR}" == "xTRUE" ]]; then

	# Define log file
	LOG="${ROOT}/logs/compile.${RUNDIR}.log"
	rm -f ${LOG}

	# Configure and build GEOS-Chem source code
	# and increment pass/fail/remain counters
	config_and_build ${ROOT} ${RUNDIR} ${LOG} ${RESULTS}
	if [[ $? -eq 0 ]]; then
	    let PASSED++
	else
	    let FAILED++
	fi
	let REMAIN--
    fi
done

#============================================================================
# Check the number of simulations that have passed
#============================================================================

# Print summary to log
print_to_log " "                                           ${RESULTS}
print_to_log "Summary of compilation test results:"        ${RESULTS}
print_to_log "${SEP_MINOR}"                                ${RESULTS}
print_to_log "Complilation tests passed:        ${PASSED}" ${RESULTS}
print_to_log "Complilation tests failed:        ${FAILED}" ${RESULTS}
print_to_log "Complilation tests not completed: ${REMAIN}" ${RESULTS}

# Check if all tests passed
if [[ "x${PASSED}" == "x${NUM_TESTS}" ]]; then
    print_to_log ""                                        ${RESULTS}
    print_to_log "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" ${RESULTS}
    print_to_log "%%%  All compilation tests passed!  %%%" ${RESULTS}
    print_to_log "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" ${RESULTS}
fi

#============================================================================
# Cleanup and quit
#============================================================================

# Free local variables
unset FAILED
unset LOG
unset NUM_TESTS
unset PASSED
unset REMAIN
unset RESULTS
unset ROOT

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
