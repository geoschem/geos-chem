#!/bin/bash

#SBATCH -n 96
#SBATCH -N 2
#SBATCH -c 1
#SBATCH --exclusive
##-------------------------------------
## Use this for 1-month runs (default)
#SBATCH -t 0-8:00
##------------------------------------
## Use this for 1-year runs
##SBATCH -t 1-0:00
##------------------------------------
#SBATCH -p sapphire,huce_cascade,seas_compute,shared
#SBATCH --mem=180000
#SBATCH --mail-type=END

###############################################################################
### Sample GCHP run script for Harvard Cannon cluster (using SLURM).
###
### -n           : Requests this many cores (across all nodes)
### -N           : Requests this number of nodes
### -c 1         : Forces 1 task per CPU (needed after 2024/03/05 SLURM update)
### --mem        : Requests this amount of memory in GB
### -p           : Requests these partitions where the job can run
### -t           : Requests time for the job (days-hours:minutes)
###  --exclusive : Reserves entire nodes (i.e. to prevent backfilling jobs)
###############################################################################

# Exit if command fails
set -e

# Debug option to print all commands executed in this script
#set -x

# Define log name to include simulation start date
start_str=$(sed 's/ /_/g' cap_restart)
log=gchp.${start_str:0:13}z.log
echo "Writing output to ${log}"

# Update config files, set restart symlink, load run env, and do sanity checks
source setCommonRunSettings.sh > ${log}
source setRestartLink.sh >> ${log}
source gchp.env >> ${log}
source checkRunSettings.sh >> ${log}

# Cannon-specific setting to get around connection issues at high # cores
export OMPI_MCL_btl=openib

# Run GCHP and evenly distribute tasks across nodes
NX=$( grep NX GCHP.rc | awk '{print $2}' )
NY=$( grep NY GCHP.rc | awk '{print $2}' )
coreCount=$(( ${NX} * ${NY} ))
planeCount=$(( ${coreCount} / ${SLURM_NNODES} ))
if [[ $(( ${coreCount} % ${SLURM_NNODES} )) > 0 ]]; then
    planeCount=$(( ${planeCount} + 1 ))
fi
time srun -n ${coreCount} -N ${SLURM_NNODES} -m plane=${planeCount} --mpi=pmix ./gchp >> ${log}

# Rename mid-run checkpoint files, if any. Discard file if time corresponds
# to run start time since duplicate with initial restart file.
chkpnts=$(ls Restarts)
for chkpnt in ${chkpnts}
do
    if [[ "$chkpnt" == *"gcchem_internal_checkpoint."* ]]; then
       chkpnt_time=${chkpnt:27:13}
       if [[ "${chkpnt_time}" = "${start_str:0:13}" ]]; then
          rm ./Restarts/${chkpnt}
       else
          new_chkpnt=./Restarts/GEOSChem.Restart.${chkpnt_time}z.c${N}.nc4
          mv ./Restarts/${chkpnt} ${new_chkpnt}
       fi
    fi
done

# Rename restart file and update restart symlink if new start time ok
new_start_str=$(sed 's/ /_/g' cap_restart)
if [[ "${new_start_str}" = "${start_str}" || "${new_start_str}" = "" ]]; then
   echo "ERROR: GCHP failed to run to completion. Check the log file for more information."
   rm -f Restarts/gcchem_internal_checkpoint
   exit 1
else
    N=$(grep "CS_RES=" setCommonRunSettings.sh | cut -c 8- | xargs )    
    mv Restarts/gcchem_internal_checkpoint Restarts/GEOSChem.Restart.${new_start_str:0:13}z.c${N}.nc4
    source setRestartLink.sh
fi

exit 0
