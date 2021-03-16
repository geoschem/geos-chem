#!/bin/bash

#SBATCH -c 24
#SBATCH -N 1
#SBATCH -t 0-00:30
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
root=$(pwd -P)

# In SLURM: Load computational environment and OpenMP settings
# Otherwise (e.g. for testing) use a small number of OpenMP threads
if [[ "x${SLURM_JOBID}" == "x" ]]; then
    export OMP_NUM_THREADS=6
else
    . ~/.bashrc
    . ${root}/gcclassic_env.sh
fi

# Load common functions for tests
. ${root}/commonFunctionsForTests.sh

# Count the number of tests to be done = number of run directories
numTests=$(ls -1 "${root}/build" | wc -l)

# All integration tests will use debugging features
baseOptions="-DCMAKE_BUILD_TYPE=Debug -DRUNDIR='' -DINSTALLCOPY=${root}/exe_files"

#============================================================================
# Initialize results logfile
#============================================================================

# Results logfile name
results="${root}/logs/results.compile.log"
rm -f ${results}

# Print header to results log file
print_to_log "${SEP_MAJOR}"                                ${results}
print_to_log "GEOS-Chem Classic: Compilation Test Results" ${results}
print_to_log ""                                            ${results}
print_to_log "Using ${OMP_NUM_THREADS} OpenMP threads"     ${results}
print_to_log "Number of compilation tests: ${numTests}"    ${results}
print_to_log "${SEP_MAJOR}"                                ${results}

#============================================================================
# Configure and compile code in each GEOS_Chem run directory
#============================================================================
print_to_log " "                   ${results}
print_to_log "Compiliation tests:" ${results}
print_to_log "${SEP_MINOR}"        ${results}

# Change to the top-level build directory
cd ${root}

# Keep track of the number of tests that passed & failed
let passed=0
let failed=0
let remain=${numTests}

# Loop over build directories
for dir in default apm bpch rrtmg tomas15 tomas40; do

    # Define build directory
    buildDir="${root}/build/${dir}"

    # Define log file
    log="${root}/logs/compile.${dir}.log"
    rm -f ${log}

    # Configure and build GEOS-Chem source code
    # and increment pass/fail/remain counters
    build_gcclassic ${root} ${buildDir} ${log} ${results} "${baseOptions}"
    if [[ $? -eq 0 ]]; then
	let passed++
    else
	let failed++
    fi
    let remain--

done

#============================================================================
# Check the number of simulations that have passed
#============================================================================

# Print summary to log
print_to_log " "                                           ${results}
print_to_log "Summary of compilation test results:"        ${results}
print_to_log "${SEP_MINOR}"                                ${results}
print_to_log "Complilation tests passed:        ${passed}" ${results}
print_to_log "Complilation tests failed:        ${failed}" ${results}
print_to_log "Complilation tests not completed: ${remain}" ${results}

# Check if all tests passed
if [[ "x${passed}" == "x${numTests}" ]]; then
    print_to_log ""                                        ${results}
    print_to_log "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" ${results}
    print_to_log "%%%  All compilation tests passed!  %%%" ${results}
    print_to_log "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" ${results}
fi

#============================================================================
# Cleanup and quit
#============================================================================

# Free local variables
unset baseOptions
unset failed
unset dir
unset log
unset numTests
unset passed
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
