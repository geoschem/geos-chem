#!/usr/bin/env bash
#===============================================================================
#SBATCH --job-name=GCHP_JOB_NAME
#SBATCH --mail-user=USER@york.ac.uk
#SBATCH --mail-type=ALL
#SBATCH --account=YOUR-VIKING-ACCOUNT-CODE
#SBATCH --output=%x_%j.log
#-------------------------------------------------------------------------------
#SBATCH --partition=nodes
#SBATCH --nodes=3
#SBATCH --exclusive
#SBATCH --ntasks-per-node=40
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --time=01:00:00
#===============================================================================
set -e

# SET UP CLEAN MODULES:
source gchp.env

# MPI PERFORMANCE TWEAKS:
export I_MPI_ADJUST_GATHERV="3"
export I_MPI_ADJUST_ALLREDUCE="12"

# CHANGE USER LIMITS:
ulimit -s unlimited

# RUN GCHP:
mpiexec -n "${SLURM_NTASKS}" ./gchp
