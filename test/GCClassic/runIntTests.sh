#!/bin/bash

#SBATCH -c 24
#SBATCH -N 1
#SBATCH -t 0-3:00
#SBATCH -p huce_cascade
#SBATCH --mem=16000
#SBATCH --mail-type=END

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: runIntTests
#
# !DESCRIPTION: Script to run the various integration tests.
#  Output is sent to the "logs" folder.
#\\
#\\
# !CALLING SEQUENCE:
#  ./runIntTests
#
# !REMARKS:
#
# !REVISION HISTORY:
#  09 Oct 2020 - R. Yantosca - Initial version
#  See the subsequent Git history with the gitk browser!
#EOP
#------------------------------------------------------------------------------
#BOC

# Separator line
line="\n====================================================\n"

# Get the long path of this folder
root=`pwd -P`

# Loop over run directories
for rundir in ./*; do

    # Loop over run directories
    if [[ -d ${rundir} && "x${rundir}" != "x./logs" ]]; then

	# Create log file for each run directory in the logs/folder
	log=${root}/logs/${rundir}.log

	# Echo information
	printf "${line}Now testing ${rundir}${line}" >> ${log} 2>&1

	# Code configuration
	cd ${rundir}/build
	cmake ../CodeDir >> ${log} 2>&1
	if [[ $? -ne 0 ]]; then
	    msg="Configuration of code in ${rundir} halted with error!\n"
	    printf ${msg} >> ${log} 2>&1
	fi

	# Code compilation
	make -j install >> ${log} 2>&1
	if [[ $? -ne 0 ]]; then
	    msg="Compilation of code in ${rundir} halted with error!\n"
	    printf ${msg} >> ${log} 2>&1
	fi

	# Code execution
	cd ../
	./gcclassic >> ${log} 2>&1
	if [[ $? -ne 0 ]]; then
	    msg="Program execution in ${rundir} halted with error!\n"
	    printf ${msg} >> {$log} 2>&1
	fi

	# Navigate to integration test top-level folder
	cd ../
    fi
done
