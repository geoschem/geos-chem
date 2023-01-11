#!/usr/bin/env bash
#===============================================================================
#SBATCH --job-name=GC-CLASSIC_JOB_NAME
#SBATCH --mail-user=USER@york.ac.uk
#SBATCH --mail-type=ALL
#SBATCH --account=YOUR-VIKING-ACCOUNT-CODE
#SBATCH --output=%x_%j.log
#-------------------------------------------------------------------------------
#SBATCH --partition=nodes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=40G
#SBATCH --time=JOB_TIME
#===============================================================================
set -e

# SET UP CLEAN MODULES:
source runScriptSamples/operational_examples/york_viking/gcclassic.intel.2021b.viking.env

# SET OPENMP VARIABLES:
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OMP_STACKSIZE=500m

# CHANGE USER LIMITS:
ulimit -s unlimited

# RUN GC CLASSIC:
srun ./gcclassic
