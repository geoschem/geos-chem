#!/bin/bash

#SBATCH -n 48
#SBATCH -N 2
#SBATCH -t 0-0:10
#SBATCH -p huce_cascade
#SBATCH --mem=110000
#SBATCH --mail-type=ALL

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

# Define GEOS-Chem log file
log="gchp.log"

# Update config files, set restart symlink, load run env, and do sanity checks
source setCommonRunSettings.sh > ${log}
source setRestartLink.sh >> ${log}
source gchp.env >> ${log}
source checkRunSettings.sh >> ${log}

if [[ $? == 0 ]]; then

    # Source your environment file. This requires first setting the gchp.env
    # symbolic link using script setEnvironment in the run directory. 
    # Be sure gchp.env points to the same file for both compilation and run.
    gchp_env=$(readlink -f gchp.env)
    if [ ! -f ${gchp_env} ] 
    then
       echo "Set symbolic link to env file using ./setEnvironment.sh."
       echo "Exiting."
       exit 1
    fi
    echo " " >> ${log}
    echo "WARNING: You are using environment settings in ${gchp_env}" >> ${log}
    source ${gchp_env} >> ${log}

    # Use SLURM to distribute tasks across nodes
    NX=$( grep NX GCHP.rc | awk '{print $2}' )
    NY=$( grep NY GCHP.rc | awk '{print $2}' )
    coreCount=$(( ${NX} * ${NY} ))
    planeCount=$(( ${coreCount} / ${SLURM_NNODES} ))
    if [[ $(( ${coreCount} % ${SLURM_NNODES} )) > 0 ]]; then
	${planeCount}=$(( ${planeCount} + 1 ))
    fi

    # Echo info from computational cores to log file for displaying results
    echo "# of CPUs : ${coreCount}" >> ${log}
    echo "# of nodes: ${SLURM_NNODES}" >> ${log}
    echo "-m plane  : ${planeCount}" >> ${log}
    echo ' ' >> ${log}

    # Harvard Cannon-specific setting to get around connection issues at high # cores
    export OMPI_MCL_btl=openib

    # Start the simulation
    echo '===> Run started at' `date` >> ${log}
    time srun -n ${coreCount} -N ${SLURM_NNODES} -m plane=${planeCount} --mpi=pmix ./gchp >> ${log}
    echo '===> Run ended at' `date` >> ${log}

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

else
    cat ${log}
fi

exit 0

