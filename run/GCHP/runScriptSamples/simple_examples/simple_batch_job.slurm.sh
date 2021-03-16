#!/bin/bash
#
#SBATCH -n 60
#SBATCH -N 2
#SBATCH -t 2:00:00
#SBATCH -p huce_intel
#SBATCH --mem=110G
#SBATCH --mail-type=ALL
#
# DESCRIPTION      Sample GCHP batch job script for Cannon
# SCHEDULER        Slurm
# DATE             2021-02-19
# SIMULATION       C90, 2019-01-01 to 2019-01-02
# RESOURCES        2 nodes, 30 processes per node
# SEE ALSO         sbatch man pages, http://gchp.rtfd.io/
#

# Load modules (so software dependencies are available)
module purge
module load intel/18.0.5-fasrc01
module load openmpi/4.0.1-fasrc01
module load netcdf-fortran/4.5.2-fasrc01
module load cmake/3.16.1-fasrc01

# Misc. configuration
module list     # print loaded modules
set -e          # if a subsequent command fails, treat it as fatal (don't continue)
set -x          # for remainder of script, echo commands to the job's log file

# Unlimit resources (to prevent OS killing GCHP due to resource usage)
ulimit -c 0                  # coredumpsize
ulimit -l unlimited          # memorylocked
ulimit -u 50000              # maxproc
ulimit -v unlimited          # vmemoryuse
ulimit -s unlimited          # stacksize

# Simple GCHP launch
rm -f cap_restart gcchem*                        # delete checkpoint/restart spec/data
./runConfig.sh                                   # update configuration files
 srun -n 60 -N 2 -m plane=30 --mpi=pmix ./gchp   # launch 60 GCHP processes