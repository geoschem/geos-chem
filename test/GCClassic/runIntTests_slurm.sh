#!/bin/bash

#SBATCH -c 24
#SBATCH -N 1
#SBATCH -t 0-12:00
#SBATCH -p huce_cascade
#SBATCH --mem=16000
#SBATCH --mail-type=END

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
#  sbatch ./runIntTests_slurm.sh
#
# !REVISION HISTORY:
#  09 Oct 2020 - R. Yantosca - Initial version
#  See the subsequent Git history with the gitk browser!
#EOP
#------------------------------------------------------------------------------
#BOC

#============================================================================
# Global variable and function definitions
#============================================================================

# Set the number of OpenMP threads to what we requested in SLURM
if [[ "x$SLURM_CPUS_PER_TASK" != "x" ]]; then
    export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
fi

# Get the long path of this folder
ROOT=`pwd -P`

# Log file with results
RESULTS=${ROOT}/logs/results.log

# Load common functions
. ${ROOT}/commonFunctionsForTests.sh

# Print the number of OpenMP threads
printf "Using ${OMP_NUM_THREADS} OpenMP threads\n" >> ${RESULTS}

#============================================================================
# Configure and compile code in each run directory
#============================================================================
for RUNDIR in *; do

    # Configure and build the model (skip logs directory)
    if [[ -d ${RUNDIR} && "x${RUNDIR}" != "xlogs" ]]; then
       LOG=${ROOT}/logs/compile.${RUNDIR}.log
       config_and_build ${ROOT} ${RUNDIR} ${LOG}
    fi
done

#============================================================================
# Configure and compile code in each run directory
#============================================================================
for RUNDIR in *; do

    # Run GEOS-Chem Classic (skip logs directory)
    if [[ -d ${RUNDIR} && "x${RUNDIR}" != "xlogs" ]]; then
       LOG=${ROOT}/logs/run.${RUNDIR}.log
       run_gcclassic ${ROOT} ${RUNDIR} ${LOG}
    fi
done
