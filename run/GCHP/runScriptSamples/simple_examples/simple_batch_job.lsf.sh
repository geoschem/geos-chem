#!/bin/bash
#
# DESCRIPTION      Sample GCHP batch job script for Compute1 (Docker-based environments)
# SCHEDULER        LSF
# DATE             2021-02-19
# SIMULATION       C90, 2019-01-01 to 2019-01-02
# RESOURCES        2 nodes, 36 processes per node
# SEE ALSO         bsub man pages, http://gchp.rtfd.io/
#
#BSUB -q rvmartin
#BSUB -n 72
#BSUB -W 2:00
#BSUB -R "rusage[mem=120GB] span[ptile=36] select[mem < 2TB]"
#BSUB -a 'docker(geoschem/gchp:13.0.0-beta.1-13-g924e47f)
#BSUB -o lsf-%J.txt

# Load modules (so software dependencies are available)
. /spack/share/spack/setup-env.sh
module load openmpi-4.0.1-gcc-9-sdj47y5

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
rm -f cap_restart gcchem*           # delete checkpoint/restart spec/data
./runConfig.sh                      # update configuration files
mpiexec -n 72 ./gchp > runlog.txt   # launch 72 GCHP processes
