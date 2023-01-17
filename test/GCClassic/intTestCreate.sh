#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: intTestCreate.sh
#
# !DESCRIPTION: Creates integration test run directories in a user-specified
#  root folder, and copies a run script there.
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
cd ../../
superProjectDir=$(pwd -P)
cd ${superProjectDir}

# Directory where the run creation scripts are found
runDir=${geosChemDir}/run/GCClassic

# Load file with utility functions to setup configuration files
. ${geosChemDir}/test/shared/commonFunctionsForTests.sh

# Get the absolute path of the root folder
root=$(absolute_path ${root})

# Log file
log=${root}/logs/createIntTests.log

# Echo header
printf "${SEP_MAJOR}\n"
printf "Creating GCClassic Integration Tests\n"
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

# Remove everything in the test folder
cleanup_files ${root}

# Make the directory for the executables
printf "\nCreating new build and executable directories:\n"
echo " ... ${root}/exe_files"
mkdir -p ${root}/exe_files

# Make the build directories
if [[ ! -d ${root}/build ]]; then
    for dir in ${EXE_BUILD_LIST[@]}; do
	echo " ... ${root}/build/${dir}"
	mkdir -p ${root}/build/${dir}
    done
fi

# Copying the run scripts to the Integration Test root folder
printf "\nCopying run scripts to: ${root}\n"
cp -f ${envFile} ${root}/gcclassic_env.sh
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
# Create individual run directories: 2x25 - MERRA2 - 72L
#=============================================================================

dir="gc_2x25_CH4_merra2"
create_rundir "3\n1\n2\n1\n${root}\n${dir}\nn\n"          ${root} ${dir} ${log}

dir="gc_2x25_CO2_merra2"
create_rundir "4\n1\n2\n1\n${root}\n${dir}\nn\n"          ${root} ${dir} ${log}

dir="gc_2x25_aerosol_merra2"
create_rundir "2\n1\n2\n1\n${root}\n${dir}\nn\n"          ${root} ${dir} ${log}

dir="gc_2x25_fullchem_merra2"
create_rundir "1\n1\n1\n2\n1\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}

# DEBUG: Exit after creating a couple of rundirs
# if the 2nd argument is passed and not a null string
if [[ "x${short}" != "x" ]]; then
    cd ${testDir}
    exit 0
fi

#----------------------------------------------------------------------------
# NOTE: Disable most 2x25 fullchem rundirs to reduce execution time
# (bmy, 05 Jan 2020)
# 
#dir="gc_2x25_fullchem_aciduptake_merra2"
#create_rundir "1\n5\n1\n2\n1\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}
#
#dir="gc_2x25_fullchem_APM_merra2"
#create_rundir "1\n7\n1\n2\n1\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}
#
# NOTE: Benchmark @ 2 x 2.5 will fail due to no available mass tuning factor
#dir="gc_2x25_fullchem_benchmark_merra2"
#create_rundir "1\n1\n2\n1\n2\n1\n${root}\n${dir}\nn\n"    ${root} ${dir} ${log}
#
#dir="gc_2x25_fullchem_complexSOA_merra2"
#create_rundir "1\n3\n1\n1\n1\n1\n${root}\n${dir}\nn\n"    ${root} ${dir} ${log}
#
#dir="gc_2x25_fullchem_complexSOA_SVPOA_merra2"
#create_rundir "1\n3\n2\n1\n1\n1\n${root}\n${dir}\nn\n"    ${root} ${dir} ${log}
#
#dir="gc_2x25_fullchem_marinePOA_merra2"
#create_rundir "1\n4\n1\n2\n1\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}
#----------------------------------------------------------------------------

dir="gc_2x25_Hg_merra2"
create_rundir "5\n1\n2\n1\n${root}\n${dir}\nn\n"          ${root} ${dir} ${log}

dir="gc_2x25_POPs_BaP_merra2"
create_rundir "6\n1\n1\n2\n1\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}

dir="gc_2x25_tagCH4_merra2"
create_rundir "7\n1\n2\n1\n${root}\n${dir}\nn\n"          ${root} ${dir} ${log}

dir="gc_2x25_tagCO_merra2"
create_rundir "8\n1\n2\n1\n${root}\n${dir}\nn\n"          ${root} ${dir} ${log}

dir="gc_2x25_tagO3_merra2"
create_rundir "9\n1\n2\n1\n${root}\n${dir}\nn\n"          ${root} ${dir} ${log}

dir="gc_2x25_TransportTracers_merra2"
create_rundir "10\n1\n2\n1\n${root}\n${dir}\nn\n"         ${root} ${dir} ${log}

dir="gc_2x25_TransportTracers_merra2_LuoWd"
create_rundir "10\n1\n2\n1\n${root}\n${dir}\nn\n"         ${root} ${dir} ${log}

dir="gc_2x25_metals_merra2"
create_rundir "11\n1\n2\n1\n${root}\n${dir}\nn\n"         ${root} ${dir} ${log}

dir="gc_2x25_carboncycle_merra2"
create_rundir "12\n1\n2\n1\n${root}\n${dir}\nn\n"         ${root} ${dir} ${log}

#=============================================================================
# Create individual run directories: 2x25 - GEOSFP - 72L
#=============================================================================

#-----------------------------------------------------------------------------
# NOTE: Disable 2x25 tests w/ GEOS-FP to reduce testing time (bmy, 05 Jan 2022)
#
#dir="gc_2x25_CH4_geosfp"
#create_rundir "3\n2\n2\n1\n${root}\n${dir}\nn\n"          ${root} ${dir} ${log}
#
#dir="gc_2x25_CO2_geosfp"
#create_rundir "4\n2\n2\n1\n${root}\n${dir}\nn\n"          ${root} ${dir} ${log}
#
#dir="gc_2x25_aerosol_geosfp"
#create_rundir "2\n2\n2\n1\n${root}\n${dir}\nn\n"          ${root} ${dir} ${log}
#
#dir="gc_2x25_fullchem_geosfp"
#create_rundir "1\n1\n2\n2\n1\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}
#
#dir="gc_2x25_fullchem_aciduptake_geosfp"
#create_rundir "1\n5\n2\n2\n1\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}
#
#dir="gc_2x25_fullchem_APM_geosfp"
#create_rundir "1\n7\n2\n2\n1\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}
#
## NOTE: Benchmark @ 2 x 2.5 will fail due to no available tuning factor
##dir="gc_2x25_fullchem_benchmark_geosfp"
##create_rundir "1\n1\n2\n2\n2\n1\n${root}\n${dir}\nn\n"    ${root} ${dir} ${log}
#
#dir="gc_2x25_fullchem_complexSOA_geosfp"
#create_rundir "1\n3\n1\n2\n1\n1\n${root}\n${dir}\nn\n"    ${root} ${dir} ${log}
#
#dir="gc_2x25_fullchem_complexSOA_SVPOA_geosfp"
#create_rundir "1\n3\n2\n2\n1\n1\n${root}\n${dir}\nn\n"    ${root} ${dir} ${log}
#
#dir="gc_2x25_fullchem_marinePOA_geosfp"
#create_rundir "1\n4\n2\n2\n1\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}
#
#dir="gc_2x25_Hg_geosfp"
#create_rundir "5\n2\n2\n1\n${root}\n${dir}\nn\n"          ${root} ${dir} ${log}
#
#dir="gc_2x25_POPs_BaP.geosfp"
#create_rundir "6\n1\n2\n2\n1\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}
#
#dir="gc_2x25_tagCH4_geosfp"
#create_rundir "7\n2\n2\n1\n${root}\n${dir}\nn\n"          ${root} ${dir} ${log}
#
#dir="gc_2x25_tagCO_geosfp"
#create_rundir "8\n2\n2\n1\n${root}\n${dir}\nn\n"          ${root} ${dir} ${log}
#
#dir="gc_2x25_tagO3_geosfp"
#create_rundir "9\n2\n2\n1\n${root}\n${dir}\nn\n"          ${root} ${dir} ${log}
#
#dir="gc_2x25_TransportTracers_geosfp"
#create_rundir "10\n2\n2\n1\n${root}\n${dir}\nn\n"         ${root} ${dir} ${log}
#
# NOTE: The metals simulation runs from 2011-2013, the earlier part of
# which is out of the range of the GEOS-FP met fields.  Disable
# the metals simulation with GEOS-FP met for now (bmy, 07 Jul 2021)
#dir="gc_2x25_metals_geosfp"
#create_rundir "11\n2\n2\n1\n${root}\n${dir}\nn\n"         ${root} ${dir} ${log}
#
#dir="gc_2x25_TransportTracers_geosfp_LuoWd"
#create_rundir "10\n2\n2\n1\n${root}\n${dir}\nn\n"         ${root} ${dir} ${log}
#-----------------------------------------------------------------------------

#=============================================================================
# Create individual run directories: 4x5 - MERRA2 - 72L
#=============================================================================

dir="gc_4x5_CH4_merra2"
create_rundir "3\n1\n1\n1\n${root}\n${dir}\nn\n"          ${root} ${dir} ${log}

dir="gc_4x5_aerosol_merra2"
create_rundir "2\n1\n1\n1\n${root}\n${dir}\nn\n"          ${root} ${dir} ${log}

dir="gc_4x5_fullchem_merra2"
create_rundir "1\n1\n1\n1\n1\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}

dir="gc_4x5_fullchem_merra2_LuoWd"
create_rundir "1\n1\n1\n1\n1\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}

dir="gc_4x5_fullchem_aciduptake_merra2"
create_rundir "1\n5\n1\n1\n1\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}

dir="gc_4x5_fullchem_APM_merra2"
create_rundir "1\n7\n1\n1\n1\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}

dir="gc_4x5_fullchem_benchmark_merra2"
create_rundir "1\n2\n1\n1\n1\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}

dir="gc_4x5_fullchem_complexSOA_merra2"
create_rundir "1\n3\n1\n1\n1\n1\n${root}\n${dir}\nn\n"    ${root} ${dir} ${log}

dir="gc_4x5_fullchem_complexSOA_SVPOA_merra2"
create_rundir "1\n3\n2\n1\n1\n1\n${root}\n${dir}\nn\n"    ${root} ${dir} ${log}

dir="gc_4x5_fullchem_marinePOA_merra2"
create_rundir "1\n4\n1\n1\n1\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}

dir="gc_4x5_fullchem_RRTMG_merra2"
create_rundir "1\n8\n1\n1\n1\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}

dir="gc_4x5_fullchem_TOMAS15_merra2_47L"
create_rundir "1\n6\n1\n1\n1\n2\n${root}\n${dir}\nn\n"    ${root} ${dir} ${log}

dir="gc_4x5_fullchem_TOMAS40_merra2_47L"
create_rundir "1\n6\n2\n1\n1\n2\n${root}\n${dir}\nn\n"    ${root} ${dir} ${log}

dir="gc_4x5_Hg_merra2"
create_rundir "5\n1\n1\n1\n${root}\n${dir}\nn\n"          ${root} ${dir} ${log}

dir="gc_4x5_POPs_BaP_merra2"
create_rundir "6\n1\n1\n1\n1\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}

dir="gc_4x5_tagCH4_merra2"
create_rundir "7\n1\n1\n1\n${root}\n${dir}\nn\n"          ${root} ${dir} ${log}

dir="gc_4x5_tagCO_merra2"
create_rundir "8\n1\n1\n1\n${root}\n${dir}\nn\n"          ${root} ${dir} ${log}

dir="gc_4x5_tagO3_merra2"
create_rundir "9\n1\n1\n1\n${root}\n${dir}\nn\n"          ${root} ${dir} ${log}

dir="gc_4x5_TransportTracers_merra2"
create_rundir "10\n1\n1\n1\n${root}\n${dir}\nn\n"         ${root} ${dir} ${log}

dir="gc_4x5_TransportTracers_merra2_LuoWd"
create_rundir "10\n1\n1\n1\n${root}\n${dir}\nn\n"         ${root} ${dir} ${log}

dir="gc_4x5_metals_merra2"
create_rundir "11\n1\n1\n1\n${root}\n${dir}\nn\n"         ${root} ${dir} ${log}

dir="gc_4x5_carboncycle_merra2"
create_rundir "12\n1\n1\n1\n${root}\n${dir}\nn\n"         ${root} ${dir} ${log}

#=============================================================================
# Create individual run directories: 4x5 - GEOSFP - 72L
#=============================================================================

dir="gc_4x5_CH4_geosfp"
create_rundir "3\n2\n1\n1\n${root}\n${dir}\nn\n"          ${root} ${dir} ${log}

dir="gc_4x5_aerosol_geosfp"
create_rundir "2\n2\n1\n1\n${root}\n${dir}\nn\n"          ${root} ${dir} ${log}

dir="gc_4x5_fullchem_geosfp"
create_rundir "1\n1\n2\n1\n1\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}

#-----------------------------------------------------------------------------
# Disable most GEOS-Chem fullchem tests w/ GEOS-FP to decrease run time
# (bmy, 05 Jan 2021)
#
#dir="gc_4x5_fullchem_geosfp_LuoWd"
#create_rundir "1\n1\n2\n1\n1\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}
#
#dir="gc_4x5_fullchem_aciduptake_geosfp"
#create_rundir "1\n5\n2\n1\n1\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}
#
#dir="gc_4x5_fullchem_APM_geosfp"
#create_rundir "1\n7\n2\n1\n1\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}
#
#dir="gc_4x5_fullchem_benchmark_geosfp"
#create_rundir "1\n2\n2\n1\n1\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}
#
#dir="gc_4x5_fullchem_complexSOA_geosfp"
#create_rundir "1\n3\n1\n2\n1\n1\n${root}\n${dir}\nn\n"    ${root} ${dir} ${log}
#
#dir="gc_4x5_fullchem_complexSOA_SVPOA_geosfp"
#create_rundir "1\n3\n2\n2\n1\n1\n${root}\n${dir}\nn\n"    ${root} ${dir} ${log}
#
#dir="gc_4x5_fullchem_marinePOA_geosfp"
#create_rundir "1\n4\n2\n1\n1\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}
#
#dir="gc_4x5_fullchem_RRTMG_geosfp"
#create_rundir "1\n8\n2\n1\n1\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}
#
#dir="gc_4x5_fullchem_TOMAS15_geosfp_47L"
#create_rundir "1\n6\n1\n2\n1\n2\n${root}\n${dir}\nn\n"    ${root} ${dir} ${log}
#
#dir="gc_4x5_fullchem_TOMAS40_geosfp_47L"
#create_rundir "1\n6\n2\n2\n1\n2\n${root}\n${dir}\nn\n"    ${root} ${dir} ${log}
#-----------------------------------------------------------------------------

dir="gc_4x5_Hg_geosfp"
create_rundir "5\n2\n1\n1\n${root}\n${dir}\nn\n"          ${root} ${dir} ${log}

dir="gc_4x5_POPs_BaP_geosfp"
create_rundir "6\n1\n2\n1\n1\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}

dir="gc_4x5_tagCH4_geosfp"
create_rundir "7\n2\n1\n1\n${root}\n${dir}\nn\n"          ${root} ${dir} ${log}

dir="gc_4x5_tagCO_geosfp"
create_rundir "8\n2\n1\n1\n${root}\n${dir}\nn\n"          ${root} ${dir} ${log}

dir="gc_4x5_tagO3_geosfp"
create_rundir "9\n2\n1\n1\n${root}\n${dir}\nn\n"          ${root} ${dir} ${log}

dir="gc_4x5_TransportTracers_geosfp"
create_rundir "10\n2\n1\n1\n${root}\n${dir}\nn\n"         ${root} ${dir} ${log}

dir="gc_4x5_TransportTracers_geosfp_LuoWd"
create_rundir "10\n2\n1\n1\n${root}\n${dir}\nn\n"         ${root} ${dir} ${log}

# NOTE: The metals simulation runs from 2011-2013, the earlier part of
# which is out of the range of the GEOS-FP met fields.  Disable
# the metals simulation with GEOS-FP met for now (bmy, 07 Jul 2021)
#dir="gc_4x5_metals_geosfp"
#create_rundir "11\n2\n1\n1\n${root}\n${dir}\nn\n"         ${root} ${dir} ${log}

dir="gc_4x5_carboncycle_geosfp"
create_rundir "12\n2\n1\n1\n${root}\n${dir}\nn\n"         ${root} ${dir} ${log}

#=============================================================================
# Create individual run directories: 4x5 and 47L (both MERRA2 and GEOSFP)
#=============================================================================

dir="gc_4x5_fullchem_merra2_47L"
create_rundir "1\n1\n1\n1\n2\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}

#-----------------------------------------------------------------------------
# Disable most GEOS-Chem fullchem_47L tests w/ GEOS-FP to decrease run time
# (bmy, 05 Jan 2021)
#dir="gc_4x5_fullchem_geosfp_47L"
#create_rundir "1\n1\n2\n1\n2\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}
#-----------------------------------------------------------------------------

#=============================================================================
# Nested-grid simulations
#=============================================================================

dir="gc_05x0625_CH4_merra2_47L_na"
create_rundir "3\n1\n3\n4\n2\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}

dir="gc_025x03125_CH4_geosfp_47L_na"
create_rundir "3\n2\n4\n4\n2\n${root}\n${dir}\nn\n"       ${root} ${dir} ${log}

dir="gc_05x0625_fullchem_merra2_47L_na"
create_rundir "1\n1\n1\n3\n4\n2\n${root}\n${dir}\nn\n"    ${root} ${dir} ${log}

dir="gc_025x03125_fullchem_geosfp_47L_na"
create_rundir "1\n1\n2\n4\n4\n2\n${root}\n${dir}\nn\n"    ${root} ${dir} ${log}

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
unset CMP_PASS_STR
unset CMP_FAIL_STR
unset EXE_PASS_STR
unset EXE_FAIL_STR
#EOC
