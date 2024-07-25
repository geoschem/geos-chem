#!/bin/bash

#SBATCH -c 8
#SBATCH -N 1
#SBATCH -t 0-0:30
#SBATCH -p REQUESTED_PARTITION
#SBATCH --mem=8000
#SBATCH --mail-type=END
#BSUB -q REQUESTED_PARTITION
#BSUB -n 8
#BSUB -W 0:30
#BSUB -R "rusage[mem=8GB] span[ptile=1] select[mem < 1TB]"
#BSUB -a 'docker(registry.gsc.wustl.edu/sleong/esm:intel-2021.1.2)'
#BSUB -o lsf-%J.txt

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: parallelTestCompile.sh
#
# !DESCRIPTION: Runs compilation tests on various GEOS-Chem Classic
#  run directories (interactively or using a scheduler).
#\\
#\\
# !CALLING SEQUENCE:
#  ./parallelTestCompile.sh        # Interactive command-line execution
#  bsub parallelTestCompile.sh     # Execution via LSF
#  sbatch parallelTestCompile.sh   # Execution via SLURM
#EOP
#------------------------------------------------------------------------------
#BOC

#============================================================================
# Global variable and function definitions
#============================================================================

# This script starts executing 1 level lower than $itRoot
ptRoot=$(cd ..; pwd)

# Include global variables & functions
. "${ptRoot}/scripts/commonFunctionsForTests.sh"

# Create local convenience variables
binDir="${ptRoot}/${BIN_DIR}"
buildDir="${ptRoot}/${BUILD_DIR}"
envDir="${ptRoot}/${ENV_DIR}"
codeDir="${ptRoot}/CodeDir"
logsDir="${ptRoot}/${LOGS_DIR}"
scriptsDir="${ptRoot}/${SCRIPTS_DIR}"
site=$(get_site_name)

# Load the user-environment and the software environment
. ~/.bashrc > /dev/null 2>&1
[[ "X${site}" == "XCANNON" ]] && . ${envDir}/gcclassic.env > /dev/null 2>&1

# All parallelization tests will use debugging features
baseOptions="-DCMAKE_BUILD_TYPE=Debug -DRUNDIR='' -DINSTALLCOPY=${binDir}"

# Get the Git commit of the superproject and submodules
head_gcc=$(export GIT_DISCOVERY_ACROSS_FILESYSTEM=1; \
           git -C "${codeDir}" log --oneline --no-decorate -1)
head_gc=$(export GIT_DISCOVERY_ACROSS_FILESYSTEM=1; \
          git -C "${codeDir}/src/GEOS-Chem" log --oneline --no-decorate -1)
head_hco=$(export GIT_DISCOVERY_ACROSS_FILESYSTEM=1; \
           git -C "${codeDir}/src/HEMCO" log --oneline --no-decorate -1)

# Scheduler-specific settings
if [[ "X${SLURM_JOBID}" != "X" ]]; then

    #-----------------------
    # SLURM settings
    #-----------------------

    # Set OMP_NUM_THREADS to the same # of cores requested with #SBATCH -c
    export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

elif [[ "X${LSB_JOBID}" != "X" ]]; then

    #-----------------------
    # LSF settings
    #-----------------------

    # Set OMP_NUM_THREADS to the same # of cores requested with #BSUB -n
    export OMP_NUM_THREADS=${$LSB_DJOB_NUMPROC}

else

    #-----------------------
    # Interactive settings
    #-----------------------

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

# Count the number of tests to be done
parallelTests=("${EXE_GCC_BUILD_LIST[@]}")
numTests=${#parallelTests[@]}

#============================================================================
# Initialize results logfile
#============================================================================

# Results logfile name
results="${logsDir}/results.compile.log"
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
print_to_log ""                                            "${results}"
if [[ "X${SLURM_JOBID}" != "X" ]]; then
    print_to_log "Submitted as SLURM job: ${SLURM_JOBID}"  "${results}"
elif  [[ "X${LSB_JOBID}" != "X" ]]; then
    print_to_log "Submitted as LSF job: ${LSB_JOBID}"      "${results}"
else
    print_to_log "Submitted as interactive job"            "${results}"
fi
print_to_log "${SEP_MAJOR}"                                "${results}"

#============================================================================
# Configure and compile code in each GEOS_Chem run directory
#============================================================================
print_to_log " "                   "${results}"
print_to_log "Compiliation tests:" "${results}"
print_to_log "${SEP_MINOR}"        "${results}"

# Change to the top-level build directory
cd "${ptRoot}"

# Keep track of the number of tests that passed & failed
let passed=0
let failed=0
let remain=${numTests}

# Loop over build directories
for dir in ${parallelTests[@]}; do

    # Define build directory
    thisBuildDir="${buildDir}/${dir}"

    # Define log file
    log="${logsDir}/compile.${dir}.log"
    rm -f "${log}"

    # Configure and build GEOS-Chem Classic source code
    # and increment pass/fail/remain counters
    build_model "gcclassic"      "${ptRoot}" "${thisBuildDir}" \
                "${baseOptions}" "${log}"    "${results}"
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

# Check for success
if [[ "x${passed}" == "x${numTests}" ]]; then

    #--------------------------
    # Successful compilation
    #--------------------------
    print_to_log ""                                        "${results}"
    print_to_log "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" "${results}"
    print_to_log "%%%  All compilation tests passed!  %%%" "${results}"
    print_to_log "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" "${results}"

    # Start the interactive execution test script upon successful finish
    if [[ "X${SLURM_JOBID}" == "X" && "x${LSB_JOBID}" == "X" ]]; then
        echo ""
        echo "Compilation tests finished!"
        ${scriptsDir}/parallelTestExecute.sh &
    fi

else

    #---------------------------
    # Unsuccessful compilation
    #---------------------------
    if [[ "X${SLURM_JOBID}" == "X" && "x${LSB_JOBID}" == "X" ]]; then
       echo ""
       echo "Compilation tests failed!  Exiting..."
    fi
fi

#============================================================================
# Cleanup and quit
#============================================================================

# Free local variables
unset baseOptions
unset binDir
unset buildDir
unset codeDir
unset failed
unset dir
unset envDir
unset head_gcc
unset head_gc
unset head_hco
unset kernel
unset log
unset numTests
unset passed
unset ptRoot
unset remain
unset results
unset scriptsDir
unset scheduler

# Free imported variables
unset FILL
unset SEP_MAJOR
unset SEP_MINOR
unset CMP_PASS_STR
unset CMP_FAIL_STR
unset EXE_PASS_STR
unset EXE_FAIL_STR
#EOC
