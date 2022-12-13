#!/bin/bash

#SBATCH -c 8
#SBATCH -N 1
#SBATCH -t 0-00:30
#SBATCH -p REQUESTED_PARTITION
#SBATCH --mem=8000
#SBATCH --mail-type=END

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: intTestCompile.sh
#
# !DESCRIPTION: Runs compilation tests on various GEOS-Chem Classic
#  run directories (interactively or using a scheduler).
#\\
#\\
# !CALLING SEQUENCE:
#  ./intTestCompile.sh        # Interactive command-line execution
#  bsub intTestCompile.sh     # Execution via LSF
#  sbatch intTestCompile.sh   # Execution via SLURM
#EOP
#------------------------------------------------------------------------------
#BOC

#============================================================================
# Global variable and function definitions
#============================================================================

# Get the long path of the integration test folder
root=$(pwd -P)

# Load the user-environment and the software environment
. ~/.bashrc                > /dev/null 2>&1
. ${root}/gcclassic_env.sh > /dev/null 2>&1

# All integration tests will use debugging features
baseOptions="-DCMAKE_BUILD_TYPE=Debug -DRUNDIR='' -DINSTALLCOPY=${root}/exe_files"
# Get the Git commit of the superproject and submodules
head_gcc=$(export GIT_DISCOVERY_ACROSS_FILESYSTEM=1; \
	   git -C "./CodeDir" log --oneline --no-decorate -1)
head_gc=$(export GIT_DISCOVERY_ACROSS_FILESYSTEM=1; \
	  git -C "./CodeDir/src/GEOS-Chem" log --oneline --no-decorate -1)
head_hco=$(export GIT_DISCOVERY_ACROSS_FILESYSTEM=1; \
	   git -C "./CodeDir/src/HEMCO" log --oneline --no-decorate -1)

if [[ "x${SLURM_JOBID}" != "x" ]]; then

    #--------------------
    # Run via SLURM
    #--------------------

    # Set number of cores to those requested with #SBATCH -c
    export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

elif [[ "x${LSB_JOBID}" != "x" ]]; then

    #--------------------
    # Run via LSF (TODO)
    #--------------------
    echo "LSF support is coming soon!"
    exit 1

else

    #--------------------
    # Run interactively
    #--------------------

    # For AWS, set $OMP_NUM_THREADS to the available cores
    kernel=$(uname -r)
    [[ "x${kernel}" == "xaws" ]] && export OMP_NUM_THREADS=$(nproc)

fi

# Sanity check: Set OMP_NUM_THREADS to 6 if it is not set
# (this may happen when running interactively)
[[ "x${OMP_NUM_THREADS}" == "x" ]] && export OMP_NUM_THREADS=6

# Sanity check: Max out the OMP_STACKSIZE if it is not set
[[ "x${OMP_STACKSIZE}" == "x" ]] && export OMP_STACKSIZE=500m

#============================================================================
# Load common variables and functions for tesets
#============================================================================

# Include global variables & functions
. ${root}/commonFunctionsForTests.sh

# Count the number of tests to be done
numTests=${#EXE_BUILD_LIST[@]}

#============================================================================
# Initialize results logfile
#============================================================================

# Results logfile name
results="${root}/logs/results.compile.log"
rm -f ${results}

# Print header to results log file
print_to_log "${SEP_MAJOR}"                                "${results}"
print_to_log "GEOS-Chem Classic: Compilation Test Results" "${results}"
print_to_log ""                                            "${results}"
print_to_log "GCClassic #${head_gcc}"                      "${results}"
print_to_log "GEOS-Chem #${head_gc}"                       "${results}"
print_to_log "HEMCO     #${head_hco}"                      "${results}"
print_to_log ""                                            "${results}"
print_to_log "Using ${OMP_NUM_THREADS} OpenMP threads"     "${results}"
print_to_log "Number of compilation tests: ${numTests}"    "${results}"
print_to_log "${SEP_MAJOR}"                                "${results}"

#============================================================================
# Configure and compile code in each GEOS_Chem run directory
#============================================================================
print_to_log " "                   "${results}"
print_to_log "Compiliation tests:" "${results}"
print_to_log "${SEP_MINOR}"        "${results}"

# Change to the top-level build directory
cd ${root}

# Keep track of the number of tests that passed & failed
let passed=0
let failed=0
let remain=${numTests}

# Loop over build directories
for dir in ${EXE_BUILD_LIST[@]}; do

    # Define build directory
    buildDir="${root}/build/${dir}"

    # Define log file
    log="${root}/logs/compile.${dir}.log"
    rm -f ${log}

    # Configure and build GEOS-Chem source code
    # and increment pass/fail/remain counters
    build_gcclassic "${root}" "${buildDir}" "${log}" "${results}" "${baseOptions}"
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
print_to_log " "                                           "${results}"
print_to_log "Summary of compilation test results:"        "${results}"
print_to_log "${SEP_MINOR}"                                "${results}"
print_to_log "Complilation tests passed:        ${passed}" "${results}"
print_to_log "Complilation tests failed:        ${failed}" "${results}"
print_to_log "Complilation tests not completed: ${remain}" "${results}"

# Check if all tests passed
if [[ "x${passed}" == "x${numTests}" ]]; then
    print_to_log ""                                        "${results}"
    print_to_log "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" "${results}"
    print_to_log "%%%  All compilation tests passed!  %%%" "${results}"
    print_to_log "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" "${results}"
fi

#============================================================================
# Cleanup and quit
#============================================================================

# Free local variables
unset baseOptions
unset failed
unset dir
unset head_gcc
unset head_gc
unset head_hco
unset kernel
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
unset CMP_PASS_STR
unset CMP_FAIL_STR
unset EXE_PASS_STR
unset EXE_FAIL_STR
#EOC
