#!/bin/bash

#BSUB -q rvmartin
#BSUB -n 24
#BSUB -W 168:00
#BSUB -R "rusage[mem=300GB] span[ptile=36] select[mem < 2TB]"
#BSUB -a 'docker(registry.gsc.wustl.edu/sleong/esm:intel-2021.1.2)'
#BSUB -o lsf-%J.txt
#
#
# source bashrc
      . /etc/bashrc
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
# LAUNCH GCHP 
#
#
    mpiexec -n 24 ./gchp > ${log}
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
