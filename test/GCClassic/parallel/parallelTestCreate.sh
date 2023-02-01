#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: paralleltestCreate.sh
#
# !DESCRIPTION: Creates parallelization test run directories in a user-specified
#  root folder, and copies a run script there.
#\\
#\\
# !CALLING SEQUENCE:
#  ./paralleltestCreate.sh /path/to/int/test/root /path/to/env-file
#  ./paralleltestCreate.sh /path/to/int/test/root /path/to/env-file quick=1
#
# !REVISION HISTORY:
#  See the subsequent Git history with the gitk browser!
#EOP
#------------------------------------------------------------------------------
#BOC

#=============================================================================
# Arguments
#=============================================================================

# Parallelization test root folder
ptRoot="${1}"
if [[ "x${ptRoot}" == "x" ]]; then
    echo "ERROR: The root-level directory for tests has not been specified!"
    exit 1
fi

# Environment file
envFile="${2}"
if [[ "x${envFile}" == "x" ]]; then
    echo "ERROR: The enviroment file (w/ module loads) has not been specified!"
    exit 1
fi
if [[ ! -f ${envFile} ]]; then
    echo "ERROR: The enviroment file is not a valid file!"
    exit 1
fi

# Run a short parallelization test?
quick="${3}"

#=============================================================================
# Global variable and function definitions
#=============================================================================

# Current directory
thisDir=$(pwd -P)
cd "${thisDir}"

# GCClassic superproject directory
cd ../../../../../
superProjectDir=$(pwd -P)
cd ${superProjectDir}

# GEOS-Chem and HEMCO submodule directories
geosChemDir="${superProjectDir}/src/GEOS-Chem"
hemcoDir="${superProjectDir}/src/HEMCO"

# Get the Git commit of the superproject and submodules
head_gcc=$(export GIT_DISCOVERY_ACROSS_FILESYSTEM=1; \
	   git -C "${superProjectDir}" log --oneline --no-decorate -1)
head_gc=$(export GIT_DISCOVERY_ACROSS_FILESYSTEM=1; \
	  git -C "${geosChemDir}" log --oneline --no-decorate -1)
head_hco=$(export GIT_DISCOVERY_ACROSS_FILESYSTEM=1; \
	   git -C "${hemcoDir}" log --oneline --no-decorate -1)

# Source the script containing utility functions and variables
. "${geosChemDir}/test/shared/commonFunctionsForTests.sh"

# Echo header
printf "${SEP_MAJOR}\n"
printf "Creating GEOS-Chem Classic Parallelization Tests\n\n"
printf "GCClassic #${head_gcc}\n"
printf "GEOS_Chem #${head_gc}\n"
printf "HEMCO     #${head_hco}\n"
printf "${SEP_MAJOR}\n"

#=============================================================================
# Initial setup of parallelization test directory
#=============================================================================

# Create parallelization test root folder if it doesn't exist
ptRoot=$(absolute_path "${ptRoot}")
[[ ! -d "${ptRoot}" ]] && mkdir -p "${ptRoot}"

# Remove everything in the parallelization test root folder
cleanup_files "${ptRoot}"

# Make the directory for the executables
printf "\nCreating new build and executable directories:\n"
echo " ... ${ptRoot}/exe_files"
mkdir -p "${ptRoot}/exe_files"

# Make the build directories
if [[ ! -d "${ptRoot}/build" ]]; then
    for dir in ${EXE_GCC_BUILD_LIST[@]}; do
	echo " ... ${ptRoot}/build/${dir}"
	mkdir -p "${ptRoot}/build/${dir}"
    done
fi

# Copying the run scripts to the parallelization test root folder
printf "\nCopying run scripts to: ${ptRoot}\n"
cp -f ${envFile}                            ${ptRoot}/gcclassic.env
cp -f ${thisDir}/parallelTest*.sh           ${ptRoot}
cp -f ${thisDir}/commonFunctionsForTests.sh ${ptRoot}

# Create a symbolic link to the code from the Parallelization Test ptRoot folder
[[ -L ${itRoot}/CodeDir ]] && ${itRoot}/CodeDir
ln -s "${superProjectDir}" ${ptRoot}/CodeDir

# Create log directory
if [[ !(-d "${ptRoot}/logs") ]]; then
    printf "\nCreating log directory: ${ptRoot}/logs\n"
    mkdir "${ptRoot}/logs"
fi

# Log file with echoback from rundir creation
log="${ptRoot}/logs/createParallelTests.log"

# Switch to folder where rundir creation scripts live
cd "${geosChemDir}/run/GCClassic"

# Change to the directory where we will create the rundirs
printf "\nCreating new run directories:\n"

#=============================================================================
# Create individual run directories: 4x5 - MERRA2 - 72L
#=============================================================================

# 4x5 merra2 CH4
create_rundir "3\n1\n1\n1\n${ptRoot}\n\nn\n" "${log}" "${ptRoot}"

# 4x5 merra2 aerosol
create_rundir "2\n1\n1\n1\n${ptRoot}\n\nn\n" "${log}" "${ptRoot}"

# DEBUG: Exit after creating a couple of rundirs if $quick is "yes"
if [[ "x${quick}" == "xyes" ]]; then
    cd ${thisDir}
    exit 0
fi

# 4x5 merra2 fullchem
create_rundir "1\n1\n1\n1\n1\n${ptRoot}\n\nn\n" "${log}" "${ptRoot}"

# 4x5 merra2 fullchem_LuoWd
dir="gc_4x5_merra2_fullchem_LuoWd"
create_rundir "1\n1\n1\n1\n1\n${ptRoot}\n${dir}\nn\n" "${log}" "${ptRoot}"

# 4x5 merra2 fullchem_aciduptake
create_rundir "1\n5\n1\n1\n1\n${ptRoot}\n\nn\n" "${log}" "${ptRoot}"

# 4x5 merra2 fullchem_APM
create_rundir "1\n7\n1\n1\n1\n${ptRoot}\n\nn\n" "${log}" "${ptRoot}"

# 4x5 merra2 fullchem_benchmark
create_rundir "1\n2\n1\n1\n1\n${ptRoot}\n\nn\n" "${log}" "${ptRoot}"

# 4x5 merra2 fullchem_complexSOA
create_rundir "1\n3\n1\n1\n1\n1\n${ptRoot}\n\nn\n" "${log}" "${ptRoot}"

# 4x5 merra2 fullchem_complexSOA_SVPOA
create_rundir "1\n3\n2\n1\n1\n1\n${ptRoot}\n\nn\n" "${log}" "${ptRoot}"

# 4x5 merra2 fullchem_marinePOA
create_rundir "1\n4\n1\n1\n1\n${ptRoot}\n\nn\n" "${log}" "${ptRoot}"

# 4x5 merra2 fullchem_RRTMG
create_rundir "1\n8\n1\n1\n1\n${ptRoot}\n\nn\n" "${log}" "${ptRoot}"

# 4x5 merra2 fullchem_TOMAS15_47L
create_rundir "1\n6\n1\n1\n1\n2\n${ptRoot}\n\nn\n" "${log}" "${ptRoot}"

# 4x5 merra2 Hg
create_rundir "5\n1\n1\n1\n${ptRoot}\n\nn\n" "${log}" "${ptRoot}"

# 4x5 merra2 POPs_BaP
create_rundir "6\n1\n1\n1\n1\n${ptRoot}\n\nn\n" "${log}" "${ptRoot}"

# 4x5 merra2 tagCH4
create_rundir "7\n1\n1\n1\n${ptRoot}\n\nn\n" "${log}" "${ptRoot}"

# 4x5 merra2 tagCO
create_rundir "8\n1\n1\n1\n${ptRoot}\n\nn\n" "${log}" "${ptRoot}"

# 4x5 merra2 tagO3
create_rundir "9\n1\n1\n1\n${ptRoot}\n\nn\n" "${log}" "${ptRoot}"

# 4x5 merra2 TransportTracers
create_rundir "10\n1\n1\n1\n${ptRoot}\n\nn\n" "${log}" "${ptRoot}"

# 4x5 merra2 TransportTracers_LuoWd
dir="gc_4x5_merra2_TransportTracers_LuoWd"
create_rundir "10\n1\n1\n1\n${ptRoot}\n${dir}\nn\n" "${log}" "${ptRoot}"

# 4x5 merra2 metals"
create_rundir "11\n1\n1\n1\n${ptRoot}\n\nn\n" "${log}" "${ptRoot}"

##4x5 merra2 carboncycle -- COMING SOON
#create_rundir "12\n1\n1\n1\n${ptRoot}\n\nn\n" "${log}" "${ptRoot}"

#=============================================================================
# Create individual run directories: 4x5 and 47L (MERRA2)
#=============================================================================

# 4x5 merra2 fullchem_47L"
create_rundir "1\n1\n1\n1\n2\n${ptRoot}\n\nn\n" "${log}" "${ptRoot}"

#=============================================================================
# Nested-grid simulations
#=============================================================================

# 05x0625 merra2 CH4_47L_na"
create_rundir "3\n1\n3\n4\n2\n${ptRoot}\n\nn\n" "${log}" "${ptRoot}"

#=============================================================================
# Cleanup and quit
#=============================================================================

# Switch back to the present directory
cd "${thisDir}"

# Free local variables
unset dir
unset geosChemDir
unset log
unset ptRoot
unset superProjectDir
unset thisDir

# Free imported variables
unset FILL
unset SEP_MAJOR
unset SEP_MINOR
unset CMP_PASS_STR
unset CMP_FAIL_STR
unset EXE_PASS_STR
unset EXE_FAIL_STR
#EOC
