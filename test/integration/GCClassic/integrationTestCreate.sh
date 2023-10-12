#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: integrationTestCreate.sh
#
# !DESCRIPTION: Creates integration test run directories in a user-specified
#  root folder, and copies a run script there.
#\\
#\\
# !CALLING SEQUENCE:
#  ./integrationTestCreate.sh /path/to/int/test/root /path/to/env-file
#  ./integrationTestCreate.sh /path/to/int/test/root /path/to/env-file quick=1
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
if [[ ! -f ${envFile} ]]; then
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
cd "${thisDir}"

# GCClassic superproject directory (absolute paths)
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
commonFuncs="${geosChemDir}/test/shared/commonFunctionsForTests.sh"
. "${commonFuncs}"

# Echo header
printf "${SEP_MAJOR}\n"
printf "Creating GEOS-Chem Classic Integration Tests\n\n"
printf "GCClassic #${head_gcc}\n"
printf "GEOS-Chem #${head_gc}\n"
printf "HEMCO     #${head_hco}\n"
printf "${SEP_MAJOR}\n"

#=============================================================================
# Create integration test folder and subdirectories
#=============================================================================

# Create integration test root folder if it doesn't exist
itRoot=$(absolute_path "${itRoot}")
[[ ! -d "${itRoot}" ]] && mkdir -p "${itRoot}"

# Create local convenience variables
binDir="${itRoot}/${BIN_DIR}"
buildDir="${itRoot}/${BUILD_DIR}"
envDir="${itRoot}/${ENV_DIR}"
execDir="${itRoot}/${EXEC_DIR}"
logsDir="${itRoot}/${LOGS_DIR}"
scriptsDir="${itRoot}/${SCRIPTS_DIR}"
rundirsDir="${itRoot}/${RUNDIRS_DIR}"

# Get absolute path of the environment file
envFile=$(absolute_path "${envFile}")

# Remove run directories in the test folder
cleanup_files "${itRoot}"

# Subdir for CMake builds (note: will create ${itRoot}
printf "\nCreating CMake build directories:\n"
for dir in ${EXE_GCC_BUILD_LIST[@]}; do
    printf " ... ${buildDir}/${dir}\n"
    mkdir -p "${buildDir}/${dir}"
done

# Subdir for executables
printf "\nCreating exe files directory ${binDir}\n"
mkdir -p "${binDir}"

# Subdir for env files
printf "Creating env files directory ${envDir}\n"
mkdir -p "${envDir}"

# Subdir for log files
printf "Creating logs directory      ${logsDir}\n"
mkdir -p "${logsDir}"

# Subdir for scripts
printf "Creating scripts directory   ${scriptsDir}\n"
mkdir -p "${scriptsDir}"

# Subdir for run directories
printf "Creating rundirs directory   ${rundirsDir}\n"
mkdir -p "${rundirsDir}"

# Create a symbolic link to the code from the Integration Test root folder
printf "Linking to superproject      ${itRoot}/CodeDir\n"
ln -s "${superProjectDir}" ${itRoot}/CodeDir

#=============================================================================
# Copy files to the proper folders
#=============================================================================

printf "\nCopying run scripts to: ${itRoot}/${SCRIPTS_DIR}\n"
cp -f ${envFile}                     ${envDir}/gcclassic.env
cp -f ${thisDir}/integration*.sh     ${scriptsDir}
cp -f ${commonFuncs}                 ${scriptsDir}
cp -f ${thisDir}/README.md           ${scriptsDir}
cp -f ${thisDir}/README.testroot.md  ${itRoot}/README.md

# This is necessary on Compute1 to make all scripts executable
chmod 755 -R ${scriptsDir}

# Log file with echoback from rundir creation
log="${logsDir}/createIntegrationTests.log"

# Switch to folder where rundir creation scripts live
cd "${geosChemDir}/run/GCClassic"

#=============================================================================
# Create individual run directories: 4x5 - MERRA2 - 72L
#=============================================================================
printf "\nCreating new run directories:\n"

# 4x5 merra2 CH4
create_rundir "3\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

# 4x5 merra2 CO2
create_rundir "4\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

# 4x5 merra2 aerosol
create_rundir "2\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

# 4x5 merra2 carbon
create_rundir "12\n1\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

# 4x5 merra2 carbon CH4 only
dir="gc_4x5_merra2_carbon_CH4only"
create_rundir "12\n2\n1\n1\n1\n${rundirsDir}\n${dir}\nn\n" "${log}"

# 4x5 merra2 fullchem
create_rundir "1\n1\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

# DEBUG: Exit after creating a couple of rundirsDirs if $quick is "yes"
if [[ "x${quick}" == "xyes" ]]; then
    cd ${thisDir}
    exit 0
fi

# 4x5 merra2 fullchem_LuoWd
dir="gc_4x5_merra2_fullchem_LuoWd"
create_rundir "1\n1\n1\n1\n1\n${rundirsDir}\n${dir}\nn\n" "${log}"

# 4x5 merra2 fullchem_aciduptake
create_rundir "1\n5\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

# 4x5 merra2 fullchem_APM
create_rundir "1\n7\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

# 4x5 merra2 fullchem_benchmark
create_rundir "1\n2\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

# 4x5 merra2 fullchem_complexSOA
create_rundir "1\n3\n1\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

# 4x5 merra2 fullchem_complexSOA_SVPOA
create_rundir "1\n3\n2\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

# 4x5 merra2 fullchem_marinePOA
create_rundir "1\n4\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

# 4x5 merra2 fullchem_RRTMG
create_rundir "1\n8\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

# 4x5 merra2 fullchem_TOMAS15_47L
create_rundir "1\n6\n1\n1\n1\n2\n${rundirsDir}\n\nn\n" "${log}"

# 4x5 merra2 fullchem_TOMAS40_47L
create_rundir "1\n6\n2\n1\n1\n2\n${rundirsDir}\n\nn\n" "${log}"

# 4x5 merra2 Hg
create_rundir "5\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

# 4x5 merra2 POPs_BaP
create_rundir "6\n1\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

# 4x5 merra2 tagCH4
create_rundir "7\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

# 4x5 merra2 tagCO
create_rundir "8\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

# 4x5 merra2 tagO3
create_rundir "9\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

# 4x5 merra2 TransportTracers
create_rundir "10\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

# 4x5 merra2 TransportTracers_LuoWd
dir="gc_4x5_merra2_TransportTracers_LuoWd"
create_rundir "10\n1\n1\n1\n${rundirsDir}\n${dir}\nn\n" "${log}"

# 4x5 merra2 metals
create_rundir "11\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

#=============================================================================
# Create individual run directories: 4x5 and 47L (MERRA2)
#=============================================================================

# 4x5 merra2 fullchem_47L
create_rundir "1\n1\n1\n1\n2\n${rundirsDir}\n\nn\n" "${log}"

#=============================================================================
# Nested-grid simulations
#=============================================================================

# 05x0625 merra2 CH4_47L_na
create_rundir "3\n1\n3\n4\n2\n${rundirsDir}\n\nn\n" "${log}"

# 05x0625 merra2 fullchem_47L_na
create_rundir "1\n1\n1\n3\n4\n2\n${rundirsDir}\n\nn\n" "${log}"

#=============================================================================
# Cleanup and quit
#=============================================================================

# Switch back to the present directory
cd "${thisDir}"

# Free local variables
unset binDir
unset buildDir
unset commonFuncs
unset dir
unset envDir
unset geosChemDir
unset itRoot
unset log
unset logsDir
unset rundirsDir
unset superProjectDir
unset scriptsDir
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
