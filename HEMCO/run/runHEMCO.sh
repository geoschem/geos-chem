#!/bin/bash

#SBATCH -c 8
#SBATCH -N 1
#SBATCH -t 0-12:00
#SBATCH -p huce_intel
#SBATCH --mem=15000
#SBATCH --mail-type=END

# Set the proper # of threads for OpenMP
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run GEOS_Chem.  The "time" command will return CPU and wall times.
# Stdout and stderr will be directed to the log files specified above.
time ./hemco_standalone -c HEMCO_sa_Config.rc

# Exit normally
exit 0
#EOC
