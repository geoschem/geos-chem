#!/bin/bash

#SBATCH -c 24
#SBATCH -N 1
#SBATCH -t 0-6:00
#SBATCH -p REQUESTED_PARTITION
#SBATCH --mem=90000
#SBATCH --mail-type=END
#BSUB -q REQUESTED_PARTITION
#BSUB -n 24
#BSUB -W 6:00
#BSUB -R "rusage[mem=90GB] span[ptile=1] select[mem < 2TB]"
#BSUB -a 'docker(registry.gsc.wustl.edu/sleong/esm:intel-2021.1.2)'
#BSUB -o lsf-%J.txt

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: integrationTestExecute.sh
#
# !DESCRIPTION: Runs execution tests on various GEOS-Chem Classic
#  run directories (using the SLURM scheduler).
#\\
#\\
# !CALLING SEQUENCE:
#  ./integrationTestExecute.sh        # Interactive command-line execution
#  bsub integrationTestExecute.sh     # Execution via LSF
#  sbatch integrationTestExecute.sh   # Execution via SLURM
#EOP
#------------------------------------------------------------------------------
#BOC

#============================================================================
# Global variable and function definitions
#============================================================================

# This script starts executing 1 level lower than $itRoot
itRoot=$(cd ..; pwd)

# Include global variables & functions
. "${itRoot}/scripts/commonFunctionsForTests.sh"

# Create local convenience variables
binDir="${itRoot}/${BIN_DIR}"
envDir="${itRoot}/${ENV_DIR}"
codeDir="${itRoot}/CodeDir"
logsDir="${itRoot}/${LOGS_DIR}"
rundirsDir="${itRoot}/${RUNDIRS_DIR}"
site=$(get_site_name)

# Load the environment and the software environment
. ~/.bashrc > /dev/null 2>&1
[[ "X${site}" == "XCANNON" ]] && . ${envDir}/gcclassic.env > /dev/null 2>&1

# Get the Git commit of the superproject and submodules
head_gcc=$(export GIT_DISCOVERY_ACROSS_FILESYSTEM=1; \
           git -C "${codeDir}" log --oneline --no-decorate -1)
head_gc=$(export GIT_DISCOVERY_ACROSS_FILESYSTEM=1; \
          git -C "${codeDir}/src/GEOS-Chem" log --oneline --no-decorate -1)
head_hco=$(export GIT_DISCOVERY_ACROSS_FILESYSTEM=1; \
           git -C "${codeDir}/src/HEMCO" log --oneline --no-decorate -1)

# Site-specific settings
if [[ "X${site}" == "XCANNON" && "X${SLURM_JOBID}" != "X" ]]; then

    #----------------------------------
    # SLURM settings (Harvard Cannon)
    #----------------------------------

    # Set OMP_NUM_THREADS to the same # of cores requested with #SBATCH -c
    export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

elif [[ "X${site}" == "XCOMPUTE1" && "X${LSB_JOBID}" != "X" ]]; then

    #---------------------------------
    # LSF settings (WashU Compute1)
    #---------------------------------

    # Set OMP_NUM_THREADS to the same # of cores requested with #BSUB -n
    export OMP_NUM_THREADS=${LSB_DJOB_NUMPROC}

else

    #---------------------------------
    # Interactive settings
    #---------------------------------
    echo ""
    echo "Execution tests running..."

    # For AWS, set $OMP_NUM_THREADS to the available cores
    kernel=$(uname -r)
    [[ "X${kernel}" == "Xaws" ]] && export OMP_NUM_THREADS=$(nproc)

fi

# Sanity check: Set OMP_NUM_THREADS to 8 if it is not set
# (this may happen when running interactively)
[[ "X${OMP_NUM_THREADS}" == "X" ]] && export OMP_NUM_THREADS=8

# Sanity check: Max out the OMP_STACKSIZE if it is not set
[[ "X${OMP_STACKSIZE}" == "X" ]] && export OMP_STACKSIZE=500m

# Count the number of tests to be run (same as the # of run directories)
numTests=$(count_rundirs "${rundirsDir}")

#============================================================================
# Initialize results logfile
#============================================================================

# Results logfile name
results="${logsDir}/results.execute.log"
rm -f "${results}"

# Print header to results log file
print_to_log "${SEP_MAJOR}"                                "${results}"
print_to_log "GEOS-Chem Classic: Execution Test Results"   "${results}"
print_to_log ""                                            "${results}"
print_to_log "GCClassic #${head_gcc}"                      "${results}"
print_to_log "GEOS-Chem #${head_gc}"                       "${results}"
print_to_log "HEMCO     #${head_hco}"                      "${results}"
print_to_log ""                                            "${results}"
print_to_log "Using ${OMP_NUM_THREADS} OpenMP threads"     "${results}"
print_to_log "Number of execution tests: ${numTests}"      "${results}"
print_to_log ""                                            "${results}"
if [[ "X${SLURM_JOBID}" != "X" ]]; then
    print_to_log "Submitted as SLURM job: ${SLURM_JOBID}"  "${results}"
elif  [[ "X${LSB_JOBID}" == "XCOMPUTE1" ]]; then
    print_to_log "Submitted as LSF job: ${LSB_JOBID}"      "${results}"
else
    print_to_log "Submitted as interactive job"            "${results}"
fi
print_to_log "${SEP_MAJOR}"                                "${results}"

#============================================================================
# Run the GEOS-Chem executable in each GEOS-Chem run directory
#============================================================================
print_to_log " "                 "${results}"
print_to_log "Execution tests:"  "${results}"
print_to_log "${SEP_MINOR}"      "${results}"

# Keep track of the number of tests that passed & failed
let passed=0
let failed=0
let remain=${numTests}

# Navigate to the directory containing individiual run directories
cd "${rundirsDir}"

# Loop over rundirs and run GEOS-Chem
for runDir in *; do

    # Expand rundir to absolute path
    runAbsPath="${rundirsDir}/${runDir}"

    # Do the following if for only valid GEOS-Chem run dirs
    expr=$(is_valid_rundir "${runAbsPath}")
    if [[ "X${expr}" == "XTRUE" ]]; then

        # Define log file
        log="${logsDir}/execute.${runDir}.log"
        rm -f "${log}"

        # Messages for execution pass & fail
        passMsg="$runDir${FILL:${#runDir}}.....${EXE_PASS_STR}"
        failMsg="$runDir${FILL:${#runDir}}.....${EXE_FAIL_STR}"

        # Get the executable file corresponding to this run directory
        exeFile=$(exe_name "gcclassic" "${runDir}")

        # Test if the executable exists
        if [[ -f "${binDir}/${exeFile}" ]]; then

            #----------------------------------------------------------------
            # If the executable file exists, we can do the test
            #----------------------------------------------------------------

            # Change to this run directory
            cd "${runAbsPath}"

            # Copy the executable file here
            cp -f "${binDir}/${exeFile}" .

            # Update to make sure the run directory is executable on Compute1
            if [[ "X${site}" == "XCOMPUTE1" ]]; then
                chmod 755 -R "${runAbsPath}"
            fi

            # Remove any leftover files in the run dir
            ./cleanRunDir.sh --no-interactive >> "${log}" 2>&1

            # Run the code if the executable is present.  Then update the
            # pass/fail counters and write a message to the results log file.
            if [[ "X${site}" == "XCANNON" && "X${SLURM_JOBID}" != "X" ]]; then
                srun -c ${OMP_NUM_THREADS} ./${exeFile} >> "${log}" 2>&1
            else
                ./${exeFile} >> "${log}" 2>&1
            fi

            # Determine if the job succeeded or failed
            if [[ $? -eq 0 ]]; then
                let passed++
                print_to_log "${passMsg}" "${results}"
            else
                let failed++
                print_to_log "${failMsg}" "${results}"
            fi

            # Navigate back to the folder containing run directories
            cd "${rundirsDir}"

        else

            #----------------------------------------------------------------
            # If the executable is missing, update the "fail" counter
            # and write the "failed" message to the results log file.
            #----------------------------------------------------------------
            let failed++
            if [[ "x${results}" != "x" ]]; then
                print_to_log "${failMsg}" "${results}"
            fi
        fi

        # Decrement the count of remaining tests
        let remain--
    fi
done

#============================================================================
# Check the number of simulations that have passed
#============================================================================

# Print summary to log
print_to_log " "                                            "${results}"
print_to_log "Summary of test results:"                     "${results}"
print_to_log "${SEP_MINOR}"                                 "${results}"
print_to_log "Execution tests passed: ${passed}"            "${results}"
print_to_log "Execution tests failed: ${failed}"            "${results}"
print_to_log "Execution tests not yet completed: ${remain}" "${results}"

# Check for success
if [[ "X${passed}" == "X${numTests}" ]]; then

    #--------------------------
    # Successful execution
    #--------------------------
    print_to_log ""                                         "${results}"
    print_to_log "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"    "${results}"
    print_to_log "%%%  All execution tests passed!  %%%"    "${results}"
    print_to_log "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"    "${results}"

    # Print success (if interactive)
    if [[ "X${SLURM_JOBID}" == "X" && "x${LSB_JOBID}" == "X" ]]; then
        echo ""
        echo "Execution tests finished!"
    fi

else

    #--------------------------
    # Unsuccessful execution
    #--------------------------
    if [[ "X${SLURM_JOBID}" == "X" && "X${LSB_JOBID}" == "X" ]]; then
        echo ""
        echo "Execution tests failed!  Exiting ..."
    fi
fi

#============================================================================
# Cleanup and quit
#============================================================================

# Free local variables
unset absRunPath
unset binDir
unset codeDir
unset envDir
unset exeFile
unset failed
unset failmsg
unset head_gcc
unset head_gc
unset head_hco
unset itRoot
unset log
unset logsDir
unset numTests
unset passed
unset passMsg
unset remain
unset results
unset rundirsDir
unset scheduler

# Free imported global variables
unset FILL
unset LINE
unset CMP_PASS_STR
unset CMP_FAIL_STR
unset EXE_PASS_STR
unset EXE_FAIL_STR
#EOC
