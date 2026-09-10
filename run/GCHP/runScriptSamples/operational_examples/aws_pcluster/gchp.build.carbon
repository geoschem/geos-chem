#!/bin/bash

#SBATCH -n 24
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 0-1:00
#SBATCH -p ____
#SBATCH --mem=60000
#SBATCH --job-name=build.carbon
#SBATCH --output=slurm.%x.%j.log

source gchp.env
cd build
cmake ../CodeDir -DRUNDIR=.. -DMECH=carbon
make -j
make install
cd ..

exit 0
