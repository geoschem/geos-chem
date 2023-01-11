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
#  ./intTestCreate.sh /path/to/int/test/root /path/to/env-file           or
#  ./intTestCreate.sh /path/to/int/test/root /path/to/env-file quick=1
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
itRoot="${1}"
if [[ "x${itRoot}" == "x" ]]; then
    echo "ERROR: The root-level directory for tests has not been specified!"
    exit 1
fi

# Environment file
envFile="${2}"
if [[ "x${envFile}" == "x" ]]; then
    echo "ERROR: The enviroment file (w/ module loads) has not been specified!"
    exit 1
fi
if [[ ! -f "${envFile}" ]]; then
    echo "ERROR: The enviroment file is not a valid file!"
    exit 1
fi

# Run a short integration test?
quick="${3}"

#=============================================================================
# Global variable and function definitions
#=============================================================================

# Current directory
thisDir=$(pwd -P)
cd ${thisDir}

# GCClassic superproject directory
cd ../../../../../..
superProjectDir=$(pwd -P)
cd "${superProjectDir}"

# GEOS-Chem and HEMCO submodule directories
geosChemDir="${superProjectDir}/src/GCHP_GridComp/GEOSChem_GridComp/geos-chem"
hemcoDir="${superProjectDir}/src/GCHP_GridComp/HEMCO_GridComp/HEMCO"

# Get the Git commit of the superproject and submodules
head_gchp=$(export GIT_DISCOVERY_ACROSS_FILESYSTEM=1; \
	    git -C "${superProjectDir}" log --oneline --no-decorate -1)
head_gc=$(export GIT_DISCOVERY_ACROSS_FILESYSTEM=1; \
	  git -C "${geosChemDir}" log --oneline --no-decorate -1)
head_hco=$(export GIT_DISCOVERY_ACROSS_FILESYSTEM=1; \
	   git -C "${hemcoDir}" log --oneline --no-decorate -1)

# Source the script containing utility functions and variables
. "${geosChemDir}/test/shared/commonFunctionsForTests.sh"

# Echo header
printf "${SEP_MAJOR}\n"
printf "Creating GCHP Integration Tests\n\n"
printf "GCHP      #${head_gchp}\n"
printf "GEOS_Chem #${head_gc}\n"
printf "HEMCO     #${head_hco}\n"
printf "${SEP_MAJOR}\n"

#=============================================================================
# Initial setup of integration test directory
#=============================================================================

# Create integration test root folder if it doesn't exist
itRoot=$(absolute_path "${itRoot}")
[[ ! -d "${itRoot}" ]] && mkdir -p "${itRoot}"

# Get absolute path of the environment file
envFile=$(absolute_path "${envFile}")

# Remove run directories in the te folder
cleanup_files "${itRoot}"

# Make the build directory
printf "\nCreating new build and executable directories:\n"
echo " ... ${itRoot}/exe_files"
mkdir -p "${itRoot}/exe_files"
if [[ ! -d "${itRoot}/build" ]]; then
     for dir in ${EXE_GCHP_BUILD_LIST[@]}; do
	echo " ... ${itRoot}/build/${dir}"
	mkdir -p "${itRoot}/build/${dir}"
     done
fi

# Copying the run scripts to the root test folder
printf "\nCopying run scripts to: ${itRoot}\n"
cp -f ${envFile}                            ${itRoot}/gchp.env
cp -f ${thisDir}/intTest*.sh                ${itRoot}
cp -f ${thisDir}/commonFunctionsForTests.sh ${itRoot}

# Create a symbolic link to the code from the Integration Test root folder
ln -s "${superProjectDir}" "${itRoot}/CodeDir"

# Create log directory
if [[ !(-d "${itRoot}/logs") ]]; then
    printf "\nCreating log directory: ${itRoot}/logs\n"
    mkdir "${itRoot}/logs"
fi

# Log file for containing echo-back from createRunDir.sh
log="${itRoot}/logs/createIntTests.log"

# Switch to folder where rundir creation scripts live
cd "${geosChemDir}/run/GCHP"

#=============================================================================
# Create the GCHP run directories
#=============================================================================
printf "\nCreating new run directories:\n"

# c24 geosfp TransportTracers
create_rundir "2\n1\n${itRoot}\n\nn\n" "${log}" "${itRoot}"

# c24 merra2 fullchem tagO3
create_rundir "4\n1\n${itRoot}\n\nn\n" "${log}" "${itRoot}"

# DEBUG: Exit after creating a couple of rundirs if $quick is "yes"
if [[ "x${quick}" == "xyes" ]]; then
    cd ${thisDir}
    exit 0
fi

# c24 merra2 fullchem_standard
create_rundir "1\n1\n1\n${itRoot}\n\nn\n" "${log}" "${itRoot}"

# c24 merra2 fullchem_benchmark
create_rundir "1\n2\n1\n${itRoot}\n\nn\n" "${log}" "${itRoot}"

# c24 merra2 fullchem_RRTMG
create_rundir "1\n8\n1\n${itRoot}\n\nn\n" "${log}" "${itRoot}"

#=============================================================================
# Cleanup and quit
#=============================================================================

# Switch back to the present directory
cd "${thisDir}"

# Free local variables
unset geosChemDir
unset hemcoDir
unset itRoot
unset log
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
