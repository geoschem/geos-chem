#!/bin/bash

#SBATCH -n 24
#SBATCH -N 1
#SBATCH -t 0-3:30
#SBATCH -p REQUESTED_PARTITION
#SBATCH --mem=90000
#SBATCH --mail-type=END
#BSUB -q REQUESTED_PARTITION
#BSUB -n 24
#BSUB -W 3:30
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
#  run directories (interactively or using a scheduler)
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

# Load the environment and the software environment
. ~/.bashrc          > /dev/null 2>&1
. ${envDir}/gchp.env > /dev/null 2>&1

# Get the Git commit of the superproject and submodules
head_gchp=$(export GIT_DISCOVERY_ACROSS_FILESYSTEM=1; \
           git -C "${codeDir}" log --oneline --no-decorate -1)
head_gc=$(export GIT_DISCOVERY_ACROSS_FILESYSTEM=1; \
          git -C "${codeDir}/src/GCHP_GridComp/GEOSChem_GridComp/geos-chem" \
          log --oneline --no-decorate -1)
head_hco=$(export GIT_DISCOVERY_ACROSS_FILESYSTEM=1; \
           git -C "${codeDir}/src/GCHP_GridComp/HEMCO_GridComp/HEMCO" \
           log --oneline --no-decorate -1)

# Determine the scheduler from the job ID (or lack of one)
scheduler="none"
[[ "x${SLURM_JOBID}" != "x" ]] && scheduler="SLURM"
[[ "x${LSB_JOBID}"   != "x" ]] && scheduler="LSF"

if [[ "x${scheduler}" == "xSLURM" ]]; then

    #---------------------------------
    # Run via SLURM (Harvard Cannon)
    #---------------------------------

    # Cannon-specific setting to get around connection issues at high # cores
    export OMPI_MCL_btl=openib

elif [[ "x${scheduler}" == "xLSF" ]]; then

    #---------------------------------
    # Run via LSF (WashU Compute1)
    #---------------------------------

    # Unlimit resources to prevent OS killing GCHP due to resource usage/
    # Alternatively you can put this in your environment file.
    ulimit -c 0                  # coredumpsize
    ulimit -l unlimited          # memorylocked
    ulimit -u 50000              # maxproc
    ulimit -v unlimited          # vmemoryuse
    ulimit -s unlimited          # stacksize

else

    #---------------------------------
    # Run interactively
    #---------------------------------

    # For AWS, set $OMP_NUM_THREADS to the available cores
    kernel=$(uname -r)
    [[ "x${kernel}" == "xaws" ]] && export OMP_NUM_THREADS=$(nproc)

fi

# Sanity check: Set OMP_NUM_THREADS to 6 if it is not set
# (this may happen when running interactively)
[[ "x${OMP_NUM_THREADS}" == "x" ]] && export OMP_NUM_THREADS=6

# Sanity check: Max out the OMP_STACKSIZE if it is not set
[[ "x${OMP_STACKSIZE}" == "x" ]] && export OMP_STACKSIZE=500m

# Count the number of tests to be run (same as the # of run directories)
numTests=$(count_rundirs "${rundirsDir}")

#============================================================================
# Initialize results logfile
#============================================================================

# Results logfile name
results="${logsDir}/results.execute.log"
rm -f "${results}"

# Print header to results log file
print_to_log "${SEP_MAJOR}"                               "${results}"
print_to_log "GCHP: Execution Test Results"               "${results}"
print_to_log ""                                           "${results}"
print_to_log "GCHP      #${head_gchp}"                    "${results}"
print_to_log "GEOS-Chem #${head_gc}"                      "${results}"
print_to_log "HEMCO     #${head_hco}"                     "${results}"
print_to_log ""                                           "${results}"
print_to_log "Number of execution tests: ${numTests}"     "${results}"
print_to_log ""                                           "${results}"
if [[ "x${scheduler}" == "xSLURM" ]]; then
    print_to_log "Submitted as SLURM job: ${SLURM_JOBID}" "${results}"
elif  [[ "x${scheduler}" == "xLSF" ]]; then
    print_to_log "Submitted as LSF job: ${LSB_JOBID}"     "${results}"
else
    print_to_log "Submitted as interactive job"           "${results}"
fi
print_to_log "${SEP_MAJOR}"                               "${results}"

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

    # Do the following if for only valid GCHP run dirs
    expr=$(is_gchp_rundir "${runAbsPath}")
    if [[ "x${expr}" == "xTRUE" ]]; then

        # Define log file
        log="${logsDir}/execute.${runDir}.log"
        rm -f "${log}"

        # Messages for execution pass & fail
        passMsg="$runDir${FILL:${#runDir}}.....${EXE_PASS_STR}"
        failMsg="$runDir${FILL:${#runDir}}.....${EXE_FAIL_STR}"

        # Get the executable file corresponding to this run directory
        exeFile=$(exe_name "gchp" "${runAbsPath}")

        # Test if the executable exists
        if [[ -f "${binDir}/${exeFile}" ]]; then

            #----------------------------------------------------------------
            # If the executable file exists, we can do the test
            #----------------------------------------------------------------

            # Change to the run directory
            cd "${runAbsPath}"

            # Copy the executable file here
            cp "${binDir}/${exeFile}" .

            # Update to make sure the run directory is executable
            # on Compute1.  We will later replace this test with
            # a test on the site name instead of on the scheduler.
            # TODO: Test on name rather than scheduler
            if [[ "x${scheduler}" == "xLSF" ]]; then
		chmod 755 -R "${runAbsPath}"
            fi

            # Remove any leftover files in the run dir
            ./cleanRunDir.sh --no-interactive >> "${log}" 2>&1

            # Link to the environment file
            ./setEnvironmentLink.sh "${envDir}/gchp.env"

            # Update config files, set links, load environment, sanity checks
            . setCommonRunSettings.sh >> "${log}" 2>&1
            . setRestartLink.sh       >> "${log}" 2>&1
            . gchp.env                >> "${log}" 2>&1
            . checkRunSettings.sh     >> "${log}" 2>&1

            # For safety's sake, remove restarts that weren't renamed
            for rst in Restarts; do
		if [[ "${rst}" =~ "gcchem_internal_checkpoint" ]]; then
		    rm -f "${rst}"
		fi
            done

            # Run GCHP and evenly distribute tasks across nodes
            if [[ "x${scheduler}" == "xSLURM" ]]; then

		#---------------------------------------------
		# Executing GCHP on SLURM (Harvard Cannon)
		#---------------------------------------------

		# Compute parameters for srun
		# See the gchp.run script in the folder:
		#  runScriptSamples/operational_examples/harvard_cannon
		NX=$(grep NX GCHP.rc | awk '{print $2}')
		NY=$(grep NY GCHP.rc | awk '{print $2}')
		coreCt=$(( ${NX} * ${NY} ))
		planeCt=$(( ${coreCt} / ${SLURM_NNODES} ))
		if [[ $(( ${coreCt} % ${SLURM_NNODES} )) > 0 ]]; then
		    planeCt=$(( ${planeCt} + 1 ))
		fi

		# Execute GCHP with srun
		srun -n ${coreCt} -N ${SLURM_NNODES} -m plane=${planeCt} \
		     --mpi=pmix ./${exeFile} >> "${log}" 2>&1

            elif [[ "x${scheduler}" == "xLSF" ]]; then

		#---------------------------------------------
		# Executing GCHP on LSF (WashU Compute1)
		#---------------------------------------------
		mpiexec -n 24 ./${exeFile} > "${log}" 2>&1

            else

		#---------------------------------------------
		# Executing GCHP interactively
		#---------------------------------------------
		mpirun -n 24 ./${exeFile} >> "${log}" 2>&1
            fi

            # Update pass/failed counts and write to results.log
            if [[ $? -eq 0 ]]; then

		# The run passed ...
                let passed++
                print_to_log "${passMsg}" "${results}"

		# ... so also rename the end-of-run restart file
		new_start_str=$(sed 's/ /_/g' cap_restart)
		N=$(grep "CS_RES=" setCommonRunSettings.sh | cut -c 8- | xargs )
		mv Restarts/gcchem_internal_checkpoint \
		   Restarts/GEOSChem.Restart.${new_start_str:0:13}z.c${N}.nc4
            else

		# The run failed
                let failed++
                print_to_log "${failMsg}" "${results}"

            fi

            # Change to root directory for next iteration
            cd "${rundirsDir}"

        else

            #----------------------------------------------------------------
            # If the executable is missing, update the "fail" counter
            # and write the "failed" message to the results log file.
            #----------------------------------------------------------------
            let failed++
            print_to_log "${failMsg}" "${results}"
        fi

        # Decrement the count of remaining tests
        let remain--
    fi
done

#============================================================================
# Check the number of simulations that have passed
#============================================================================

# Print summary to log
print_to_log " "                                            ${results}
print_to_log "Summary of test results:"                     ${results}
print_to_log "${SEP_MINOR}"                                 ${results}
print_to_log "Execution tests passed: ${passed}"            ${results}
print_to_log "Execution tests failed: ${failed}"            ${results}
print_to_log "Execution tests not yet completed: ${remain}" ${results}

# Check if all tests passed
if [[ "x${passed}" == "x${numTests}" ]]; then
    print_to_log ""                                         ${results}
    print_to_log "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"    ${results}
    print_to_log "%%%  All execution tests passed!  %%%"    ${results}
    print_to_log "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"    ${results}

    # Print success (if interactive)
    if [[ "x${SLURM_JOBID}" == "x" && "x${LSB_JOBID}" == "x" ]]; then
        echo ""
        echo "Execution tests finished!"
    fi

else

    #--------------------------
    # Unsuccessful execution
    #--------------------------
    if [[ "x${SLURM_JOBID}" == "x" && "x${LSB_JOBID}" == "x" ]]; then
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
unset coreCt
unset exeFile
unset failed
unset failmsg
unset head_gchp
unset head_gc
unset head_hco
unset itRoot
unset log
unset logsDir
unset numTests
unset NX
unset NY
unset passed
unset passMsg
unset planeCt
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
