#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: intTestCreate.sh
#
# !DESCRIPTION: Creates GCHP integration test run directories in a
#  user-specified root folder, and copies a run script there.
#\\
#\\
# !CALLING SEQUENCE:
#  ./intTestCreate.sh /path/to/int/test/root ENV-FILE           or
#  ./intTestCreate.sh /path/to/int/test/root ENV-FILE short=1
#
# !REMARKS:
#  Right now we pass values to the existing ./createRunDir.sh,
#  but will implement a more elegant solution later.
#
# !REVISION HISTORY:
#  09 Oct 2020 - R. Yantosca - Initial version
#  See the subsequent Git history with the gitk browser!
#EOP
#------------------------------------------------------------------------------
#BOC

#=============================================================================
# Arguments
#=============================================================================

# Integration test root folder
root=${1}
if [[ "x${root}" == "x" ]]; then
    echo "ERROR: The root-level directory for tests has not been specified!"
    exit 1
fi

# Environment file
envFile=${2}
if [[ "x${envFile}" == "x" ]]; then
    echo "ERROR: The enviroment file (w/ module loads) has not been specified!"
    exit 1
fi
if [[ ! -f ${envFile} ]]; then
    echo "ERROR: The enviroment file is not a valid file!"
    exit 1
fi

# Run a short integration test?
short=${3}

#=============================================================================
# Global variable and function definitions
#=============================================================================

# Current directory
testDir=$(pwd -P)
cd ${testDir}

# Top-level GEOS-Chem directory
cd ../..
geosChemDir=$(pwd -P)

# GCClassic superproject directory
cd ../../../..
superProjectDir=$(pwd -P)
cd ${superProjectDir}

# Directory where the run creation scripts are found
runDir=${geosChemDir}/run/GCHP

# Load file with utility functions to setup configuration files
. ${geosChemDir}/test/shared/commonFunctionsForTests.sh

# Get the absolute path of the root folder
root=$(absolute_path ${root})

# Log file
log=${root}/logs/createIntTests.log

# Echo header
printf "${SEP_MAJOR}\n"
printf "Creating GCHP Integration Tests\n"
printf "${SEP_MAJOR}\n"

#=============================================================================
# Initial setup of integration test directory
#=============================================================================

# Create integration test root folder if it doesn't exist
if [[ ! -d ${root} ]]; then
    mkdir -p ${root}
fi

# Get the absolute path of the root folder again;
# now that the folder exists
root=$(absolute_path ${root})

# Remove run directories in the test folder
cleanup_files ${root}

# Make the build directory
printf "\nCreating new build and directories:\n"
if [[ ! -d ${root}/build ]]; then
    echo " ... ${root}/build/"
    mkdir -p ${root}/build/
fi

# Copying the run scripts to the root test folder
printf "\nCopying run scripts to: ${root}\n"
cp -f ${envFile} ${root}/gchp.env
cp -f ${testDir}/intTest*.sh ${root}
cp -f ${testDir}/commonFunctionsForTests.sh ${root}

# Create a symbolic link to the code from the Integration Test root folder
ln -s ${superProjectDir} ${root}/CodeDir

# Create log directory
if [[ !(-d ${root}/logs) ]]; then
    printf "\nCreating log directory: ${root}/logs\n"
    mkdir ${root}/logs
fi

# Change to the directory where we will create the rundirs
printf "\nCreating new run directories:\n"

# Switch to folder where rundir creation scripts live
cd ${runDir}

#=============================================================================
# TransportTracers run directories
#=============================================================================

dir="gchp_TransportTracers_merra2_c24"
create_rundir "2\n1\n${root}\n${dir}\nn\n" ${root} ${dir} ${log}
ln -s ${root}/gchp.env ${root}/${dir}/gchp.env
cp ${testDir}/gchp.slurm.sh ${root}/${dir}/gchp.slurm.sh
cp ${testDir}/gchp.lsf.sh   ${root}/${dir}/gchp.lsf.sh
sed -i -e "s/CS_RES=48/CS_RES=24/" ${root}/${dir}/runConfig.sh

dir="gchp_TransportTracers_geosfp_c24"
create_rundir "2\n1\n${root}\n${dir}\nn\n" ${root} ${dir} ${log}
ln -s ${root}/gchp.env ${root}/${dir}/gchp.env
cp ${testDir}/gchp.slurm.sh ${root}/${dir}/gchp.slurm.sh
cp ${testDir}/gchp.lsf.sh   ${root}/${dir}/gchp.lsf.sh
sed -i -e "s/CS_RES=48/CS_RES=24/" ${root}/${dir}/runConfig.sh

dir="gchp_TransportTracers_merra2_c48"
create_rundir "2\n1\n${root}\n${dir}\nn\n" ${root} ${dir} ${log}
ln -s ${root}/gchp.env ${root}/${dir}/gchp.env
cp ${testDir}/gchp.slurm.sh ${root}/${dir}/gchp.slurm.sh
cp ${testDir}/gchp.lsf.sh   ${root}/${dir}/gchp.lsf.sh

dir="gchp_TransportTracers_geosfp_c48"
create_rundir "2\n1\n${root}\n${dir}\nn\n" ${root} ${dir} ${log}
ln -s ${root}/gchp.env ${root}/${dir}/gchp.env
cp ${testDir}/gchp.slurm.sh ${root}/${dir}/gchp.slurm.sh
cp ${testDir}/gchp.lsf.sh   ${root}/${dir}/gchp.lsf.sh

#=============================================================================
# Benchmark run directories
#=============================================================================

dir="gchp_fullchem_standard_merra2_c24"
create_rundir "1\n1\n1\n${root}\n${dir}\nn\n" ${root} ${dir} ${log}
ln -s ${root}/gchp.env ${root}/${dir}/gchp.env
cp ${testDir}/gchp.slurm.sh ${root}/${dir}/gchp.slurm.sh
cp ${testDir}/gchp.lsf.sh   ${root}/${dir}/gchp.lsf.sh
sed -i -e "s/CS_RES=48/CS_RES=24/" ${root}/${dir}/runConfig.sh

dir="gchp_fullchem_standard_geosfp_c24"
create_rundir "1\n1\n2\n${root}\n${dir}\nn\n" ${root} ${dir} ${log}
ln -s ${root}/gchp.env ${root}/${dir}/gchp.env
cp ${testDir}/gchp.slurm.sh ${root}/${dir}/gchp.slurm.sh
cp ${testDir}/gchp.lsf.sh   ${root}/${dir}/gchp.lsf.sh
sed -i -e "s/CS_RES=48/CS_RES=24/" ${root}/${dir}/runConfig.sh

dir="gchp_fullchem_standard_merra2_c48"
create_rundir "1\n1\n1\n${root}\n${dir}\nn\n" ${root} ${dir} ${log}
ln -s ${root}/gchp.env ${root}/${dir}/gchp.env
cp ${testDir}/gchp.slurm.sh ${root}/${dir}/gchp.slurm.sh
cp ${testDir}/gchp.lsf.sh   ${root}/${dir}/gchp.lsf.sh

dir="gchp_fullchem_standard_geosfp_c48"
create_rundir "1\n1\n2\n${root}\n${dir}\nn\n" ${root} ${dir} ${log}
ln -s ${root}/gchp.env ${root}/${dir}/gchp.env
cp ${testDir}/gchp.slurm.sh ${root}/${dir}/gchp.slurm.sh
cp ${testDir}/gchp.lsf.sh   ${root}/${dir}/gchp.lsf.sh

#=============================================================================
# Benchmark run directories
#=============================================================================

dir="gchp_fullchem_benchmark_merra2_c48"
create_rundir "1\n2\n1\n${root}\n${dir}\nn\n" ${root} ${dir} ${log}
ln -s ${root}/gchp.env ${root}/${dir}/gchp.env
cp ${testDir}/gchp.slurm.sh ${root}/${dir}/gchp.slurm.sh
cp ${testDir}/gchp.lsf.sh   ${root}/${dir}/gchp.lsf.sh

dir="gchp_fullchem_benchmark_geosfp_c48"
create_rundir "1\n2\n2\n${root}\n${dir}\nn\n" ${root} ${dir} ${log}
ln -s ${root}/gchp.env ${root}/${dir}/gchp.env
cp ${testDir}/gchp.slurm.sh ${root}/${dir}/gchp.slurm.sh
cp ${testDir}/gchp.lsf.sh   ${root}/${dir}/gchp.lsf.sh

#=============================================================================
# Cleanup and quit
#=============================================================================

# Switch back to the present directory
cd ${testDir}

# Free local variables
unset root
unset testDir
unset geosChemDir
unset superProjectDir
unset runDir
unset log
unset dir

# Free imported variables
unset FILL
unset SEP_MAJOR
unset SEP_MINOR
unset SED_INPUT_GEOS_1
unset SED_INPUT_GEOS_2
unset SED_HISTORY_RC
unset CMP_PASS_STR
unset CMP_FAIL_STR
unset EXE_PASS_STR
unset EXE_FAIL_STR
#EOC
