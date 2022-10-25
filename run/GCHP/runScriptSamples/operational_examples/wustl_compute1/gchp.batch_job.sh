#!/bin/bash

# Template file to adapt for running a GCHP batch job. Examples are
# provided for 3 job schedulers: PBS, LSF, and SLURM. Uncomment the
# code for your system and further adapt as needed.

#####################################################################
# -----------------
# Example 2: LSF
# 
#BSUB -q rvmartin
#BSUB -n 24
#BSUB -W 168:00
#BSUB -R "rusage[mem=300GB] span[ptile=36] select[mem < 2TB]"
#BSUB -a 'docker(registry.gsc.wustl.edu/sleong/esm:intel-2021.1.2)'
#BSUB -o lsf-%J.txt
#
# DESCRIPTION      Sample GCHP batch job script for Compute1 (Docker-based environments)
# SCHEDULER        LSF
# DATE             2021-02-19
# SIMULATION       C90, 2019-01-01 to 2019-01-02
# RESOURCES        2 nodes, 36 processes per node
# SEE ALSO         bsub man pages, http://gchp.rtfd.io/
#
#
# LOAD LIBRARIES - 2 options
#
# Option 1: Define environment settings in external file, link to it from
#           run directory symlink gchp.env, and source it here.
# 
#    source gchp.env
#
# Option 2: Define environment settings directly in this file. Two
#           examples are below.
# 
#    Example 1: PBS on NASA Pleiades
#
#       module purge
#       module use /nobackup/lbindle/modulefiles
#       module load gchp-env/2021.06-gnu
#
#    Example 2: LSF on WashU Compute 1 (uses docker containers)
#
      . /etc/bashrc
#      . /opt/spack/setup-env.sh
#      module load openmpi-4.0.1-gcc-9-sdj47y5
#

#################################################################
#
# ADDITIONAL PRE-RUN CONFIGURATION
#
# If a subsequent command fails, treat it as fatal (don't continue)
set -e

# For remainder of script, echo commands to the job's log file
set -x

# Unlimit resources to prevent OS killing GCHP due to resource usage/
# Alternatively you can put this in your environment file.
ulimit -c 0                  # coredumpsize
ulimit -l unlimited          # memorylocked
ulimit -u 50000              # maxproc
ulimit -v unlimited          # vmemoryuse
ulimit -s unlimited          # stacksize

#################################################################
#
# PRE-RUN COMMANDS
#
# Print loaded modules
#module list     

# Define log name to include simulation start date
start_str=$(sed 's/ /_/g' cap_restart)
log=gchp.${start_str:0:13}z.log

# Update config files, set restart symlink, and do sanity checks
source setCommonRunSettings.sh
source setRestartLink.sh
source checkRunSettings.sh

#################################################################
#
# LAUNCH GCHP - 3 examples
#
# There are several commands that may work to execute GCHP,
# and the provided examples may not be appropriate for your system.
# Check with your system administrator or online documentation about
# the appropriate way to submit jobs for your scheduler and MPI.
#
# Example 1: PBS
#
#    mpiexec -n 48 ./gchp > ${log}
#
# Example 2: LSF
#
    mpiexec -n 24 ./gchp > ${log}
#
# Example 3: SLURM
#
#    srun -n 48 -N 2 -m plane=24 --mpi=pmix ./gchp > ${log}
#

#################################################################
#
# POST-RUN COMMANDS
#
# If new start time in cap_restart is okay, rename and move restart file
# and update restart symlink
new_start_str=$(sed 's/ /_/g' cap_restart)
if [[ "${new_start_str}" = "${start_str}" || "${new_start_str}" = "" ]]; then
   echo "ERROR: cap_restart either did not change or is empty."
   exit 1
else
    N=$(grep "CS_RES=" setCommonRunSettings.sh | cut -c 8- | xargs )    
    mv gcchem_internal_checkpoint Restarts/GEOSChem.Restart.${new_start_str:0:13}z.c${N}.nc4
    source setRestartLink.sh
fi
