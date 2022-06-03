#!/bin/bash

#SBATCH -n 12
#SBATCH -N 1
#SBATCH -t 0-1:00
#SBATCH -p huce_intel
#SBATCH --mem=110000
#SBATCH --mail-type=END

# Resource request tips: 
#  (1) Use #SBATCH -n 6 to request 6 cores total across all nodes
#  (2) Use #SBATCH --exclusive to prevent other users from sharing nodes in use
#      (ability to use --exclusive may be disabled on some clusters)
#  (3) Use --mem=50G to request 50 Gigabytes per node
#  (4) Performance is enhanced by requesting entire nodes and all memory
#      even if GCHP is run on fewer cores per node than available.
#
# See SLURM documentation for descriptions of all possible settings.
# Type 'man sbatch' at the command prompt to browse documentation.

# Exit if an error is encountered
set -e

# Define GEOS-Chem log file (and remove any prior log file)
thisDir=$(pwd -P)
root=$(dirname ${thisDir})
runDir=$(basename ${thisDir})
log="${root}/logs/execute.${runDir}.log"
rm -f ${log}

# Save start date string from cap_restart
start_str=$(echo $(cat cap_restart) | sed 's/ /_/g')

# Update config files, set restart symlink, load run env, and do sanity checks
source setCommonRunSettings.sh > ${log}
source setRestartLink.sh >> ${log}
source gchp.env >> ${log}
source checkRunSettings.sh >> ${log}

# Harvard Cannon-specific setting to avoid connection issues at high # cores
export OMPI_MCL_btl=openib

# Use SLURM to distribute tasks across nodes
NX=$( grep NX GCHP.rc | awk '{print $2}' )
NY=$( grep NY GCHP.rc | awk '{print $2}' )
coreCount=$(( ${NX} * ${NY} ))
planeCount=$(( ${coreCount} / ${SLURM_NNODES} ))
if [[ $(( ${coreCount} % ${SLURM_NNODES} )) > 0 ]]; then
	${planeCount}=$(( ${planeCount} + 1 ))
fi
time srun -n ${coreCount} -N ${SLURM_NNODES} -m plane=${planeCount} --mpi=pmix ./gchp >> ${log}

# Rename and move restart file and update restart symlink if new start time ok
new_start_str=$(echo $(cat cap_restart) | sed 's/ /_/g')
if [[ "${new_start_str}" = "${start_str}" || "${new_start_str}" = "" ]]; then
   echo "ERROR: cap_restart either did not change or is empty."
   exit 1
else
    N=$(grep "CS_RES=" setCommonRunSettings.sh | cut -c 8- | xargs )    
    mv gcchem_internal_checkpoint Restarts/GEOSChem.Restart.${new_start_str}z.c${N}.nc4
    source setRestartLink.sh
fi

# Update the results log
cd ..
./intTestResults.sh
cd ${thisDir}

exit 0

