#!/bin/bash
#BSUB -q rvmartin
#BSUB -n 720
#BSUB -W 36:00
#BSUB -Q 163
#BSUB -R "rusage[mem=300GB] span[ptile=36] select[mem < 2TB]"
#BSUB -a 'docker(registry.gsc.wustl.edu/sleong/base-engineering)'
#BSUB -J "C360"
#BSUB -g /liam.bindle/C360
#BSUB -N
#BSUB -u liam.bindle@wustl.edu
#BSUB -o lsf-%J.txt

# This job auto-suspends every ~24 hours. This is to prevent blocking other jobs that are
# queued in 'rvmartin'. It relies on two files that act as semaphores:
#     CONTINUE_SEM: if exists, continue the simulation from the latest restart file
#     REQUEUE_SEM: if exists, requeue the job after it finishes
# GCHP output is logged to a file with the name <SSD>_<RTS>.log, where <SSD> is the date
# that the simulation segment is starting from, and <RTS> is the real-word start time of
# the job. 

# Source /etc/bashrc
. /etc/bashrc

# Custom functions
function last_checkpoint() {
    ls -1 gcchem_internal_checkpoint*.nc4 | tail -n 1
}
function last_checkpoint_date() {
    last_checkpoint | sed 's/gcchem_internal_checkpoint.\(20[12][0-9][0-9][0-9][0123][0-9]\).*/\1/'
}

# Set up runtime environment
set -x                           # Print executed commands
ulimit -c 0                      # coredumpsize
ulimit -l unlimited              # memorylocked
ulimit -u 50000                  # maxproc
ulimit -v unlimited              # vmemoryuse
ulimit -s unlimited              # stacksize

export TMPDIR="$__LSF_JOB_TMPDIR__"
export OMP_NUM_THREADS=1
export I_MPI_ADJUST_GATHERV=3    # IMPORTANT
export I_MPI_ADJUST_ALLREDUCE=12
export I_MPI_DAPL_UD=enable
export I_MPI_SHM_HEAP_VSIZE=512  # Might fix: Assertion failed in file ../../src/mpid/ch4/src/intel/ch4_shm_coll.c at line 2147: comm->shm_numa_layout[my_numa_node].base_addr

RESTART_DATE=SIM_START
if [ ! -f CONTINUE_SEM ] ; then  # start simulation from time=0
    rm -f cap_restart #gcchem* 
    ./setCommonRunSettings.sh
    touch CONTINUE_SEM
else                             # start simulation from the last checkpoint's date
    RESTART_FILE=$(last_checkpoint)
    RESTART_DATE=$(last_checkpoint_date)
    echo "$RESTART_DATE 000000" > cap_restart
    sed -i "s/GCHPchem_INTERNAL_RESTART_FILE: .*/GCHPchem_INTERNAL_RESTART_FILE: $RESTART_FILE/g" GCHP.rc
fi

rm -f gcchem_internal_checkpoint
mpirun -np 720 ./gchp &> ${RESTART_DATE}-$(date +"%Y%m%d_%H%M").log

if [ -f REQUEUE_SEM ] ; then
    exit 163   # exit code 163 requeues the job (according to #BSUB -Q 163)
fi
