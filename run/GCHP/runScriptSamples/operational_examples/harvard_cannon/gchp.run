#!/bin/bash

#SBATCH -n 96
#SBATCH -N 2
#SBATCH --exclusive
#SBATCH -t 0-8:00
#SBATCH -p seas_compute
#SBATCH --mem=180000
#SBATCH --mail-type=END

# Exit if command fails, and log all commands
set -xe

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
	${planeCount}=$(( ${planeCount} + 1 ))
fi
time srun -n ${coreCount} -N ${SLURM_NNODES} -m plane=${planeCount} --mpi=pmix ./gchp >> ${log}

# Rename and move restart file and update restart symlink if new start time ok
new_start_str=$(sed 's/ /_/g' cap_restart)
if [[ "${new_start_str}" = "${start_str}" || "${new_start_str}" = "" ]]; then
   echo "ERROR: GCHP failed to run to completion. Check the log file for more information."
   exit 1
else
    N=$(grep "CS_RES=" setCommonRunSettings.sh | cut -c 8- | xargs )    
    mv gcchem_internal_checkpoint Restarts/GEOSChem.Restart.${new_start_str:0:13}z.c${N}.nc4
    source setRestartLink.sh
fi

exit 0
