#!/bin/bash
#PBS -S /bin/bash
#PBS -N SGIdebugGCHP
#PBS -l select=2:ncpus=24:mpiprocs=48:model=bro
#PBS -l walltime=00:30:00
#PBS -j oe
#PBS -W group_list=[YOUR_ACCOUNT]
#PBS -m e
#PBS -M your.name@email.com

# to run on more nodes / processes, use 
# #PBS -l select=NNODEES:ncpus=NCPUS_PER_NODE:mpirocs=(NNODES * NCPUS_PER_NODE):model=bro
# for example, to run on 48 cores total over 2 nodes
# #PBS -l select=2:ncpus=24:mpiprocs=48:model=bro

# By default, PBS executes your job from your home directory.
# However, you can use the environment variable
# PBS_O_WORKDIR to change to the directory where
# you submitted your job.

cd $PBS_O_WORKDIR

# Define GEOS-Chem log file
log="gchp.log"

# Always remove cap_restart if not doing a multi-segmented run
if [[ -e cap_restart ]]; then
   rm cap_restart
fi

# Assume success until overwritten
rc=0

# Sync all config files with settings in runConfig.sh                           
./runConfig_adj.sh > ${log}
rc=$?
if [[ $rc == 0 ]]; then

    # Source your environment file. This requires first setting the gchp.env
    # symbolic link using script setEnvironment in the run directory. 
    # Be sure gchp.env points to the same file for both compilation and 
    # running. You can copy or adapt sample environment files located in 
    # ./envSamples subdirectory.
    gchp_env=$(readlink -f gchp.env)
    if [ ! -f ${gchp_env} ] 
    then
       echo "ERROR: gchp.rc symbolic link is not set!"
       echo "Copy or adapt an environment file from the ./envSamples "
       echo "subdirectory prior to running. Then set the gchp.env "
       echo "symbolic link to point to it using ./setEnvironment."
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
    
    #planeCount=$(( ${coreCount} / ${SLURM_NNODES} ))
    #if [[ $(( ${coreCount} % ${SLURM_NNODES} )) > 0 ]]; then
	#${planeCount}=$(( ${planeCount} + 1 ))
    #fi

    # Echo info from computational cores to log file for displaying results
    echo "# of CPUs : ${coreCount}" >> ${log}
    #echo "# of nodes: ${SLURM_NNODES}" >> ${log}
    echo "-m plane  : ${planeCount}" >> ${log}

    # Optionally compile 
    # Uncomment the line below to compile from scratch
    # See other compile options with 'make help'
    # make build_all

    # Echo start date
    echo ' ' >> ${log}
    echo '===> Run started at' `date` >> ${log}

    mpiexec -n $coreCount ./gchp >> $log 2>&1 &
    tail --pid=$! -f $log
    #mpiexec dplace -s1 -c 4-11 ./grinder < run_input > output
    rc=$?
# Echo end date
    echo '===> Run ended at' `date` >> ${log}
    echo "Exit code: $rc"

else
    cat ${log}
fi

# -end of script-
exit $rc
