#!/bin/bash

#SBATCH -c 24
#SBATCH -N 1
#SBATCH -t 0-03:00
#SBATCH -p REQUESTED_PARTITION
#SBATCH --mem=90000
#SBATCH --mail-type=END

#BSUB -q REQUESTED_PARTITION
#BSUB -n 24
#BSUB -W 3:00
#BSUB -R "rusage[mem=90GB] span[ptile=1] select[mem < 2TB]"
#BSUB -a 'docker(registry.gsc.wustl.edu/sleong/esm:intel-2021.1.2)'
#BSUB -o lsf-%J.txt

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
#  ./intTestExecute.sh        # Interactive command-line execution
#  bsub intTestExecute.sh     # Execution via LSF
#  sbatch intTestExecute.sh   # Execution via SLURM
#EOP
#------------------------------------------------------------------------------
#BOC

#============================================================================
# Global variable and function definitions
#============================================================================

# Get the long path of the integration test root folder
root=$(pwd -P)

# Load the environment and the software environment
. ~/.bashrc                > /dev/null 2>&1
. ${root}/gcclassic_env.sh > /dev/null 2>&1

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
    echo "LSF support coming soon!"
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

# Count the number of tests to be run (same as the # of run directories)
numTests=$(count_rundirs "${root}")

#============================================================================
# Initialize results logfile
#============================================================================

# Results logfile name
results="${root}/logs/results.execute.log"
rm -f "${results}"

# Print header to results log file
print_to_log "${SEP_MAJOR}"                               "${results}"
print_to_log "GEOS-Chem Classic: Execution Test Results"  "${results}"
print_to_log ""                                           "${results}"
print_to_log "GCClassic #${head_gcc}"                     "${results}"
print_to_log "GEOS-Chem #${head_gc}"                      "${results}"
print_to_log "HEMCO     #${head_hco}"                     "${results}"
print_to_log ""                                           "${results}"
print_to_log "Using ${OMP_NUM_THREADS} OpenMP threads"    "${results}"
print_to_log "Number of execution tests: ${numTests}"     "${results}"
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

# Loop over rundirs and run GEOS-Chem
for runDir in *; do

    # Do the following if for only valid GEOS-Chem run dirs
    expr=$(is_valid_rundir "${root}/${runDir}")
    if [[ "x${expr}" == "xTRUE" ]]; then

	# Define log file
	log="${root}/logs/execute.${runDir}.log"
	rm -f "${log}"

	# Messages for execution pass & fail
	passMsg="$runDir${FILL:${#runDir}}.....${EXE_PASS_STR}"
	failMsg="$runDir${FILL:${#runDir}}.....${EXE_FAIL_STR}"

	# Get the executable file corresponding to this run directory
	exeFile=$(gcclassic_exe_name "${runDir}")

	# Test if the executable exists
	if [[ -f "${root}/exe_files/${exeFile}" ]]; then

	    #----------------------------------------------------------------
	    # If the executable file exists, we can do the test
	    #----------------------------------------------------------------

	    # Change to this run directory; remove leftover log file
	    cd "${root}/${runDir}"

	    # Copy the executable file here
	    cp -f "${root}/exe_files/${exeFile}" .

	    # Run the code if the executable is present.  Then update the
	    # pass/fail counters and write a message to the results log file.
	    srun -c ${OMP_NUM_THREADS} "./${exeFile}" >> "${log}" 2>&1
	    if [[ $? -eq 0 ]]; then
		let passed++
		if [[ "x${results}" != "x" ]]; then
		    print_to_log "${passMsg}" "${results}"
		fi
	    else
		let failed++
		if [[ "x${results}" != "x" ]]; then
		    print_to_log "${failMsg}" "${results}"
		fi
	    fi

	    # Change to root directory for next iteration
	    cd ${root}

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
fi

#============================================================================
# Cleanup and quit
#============================================================================

# Free local variables
unset exeFile
unset failed
unset failmsg
unset head_gcc
unset head_gc
unset head_hco
unset log
unset numTests
unset passed
unset passMsg
unset remain
unset results
unset root

# Free imported global variables
unset FILL
unset LINE
unset LINELC
unset CMP_PASS_STR
unset CMP_FAIL_STR
unset EXE_PASS_STR
unset EXE_FAIL_STR
#EOC
