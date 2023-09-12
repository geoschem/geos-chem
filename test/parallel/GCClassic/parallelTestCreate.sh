#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: paralleltestCreate.sh
#
# !DESCRIPTION: Creates parallelization test run directories in a
#  user-specified root folder, and copies a run script there.
#\\
#\\
# !CALLING SEQUENCE:
#  ./paralleltestCreate.sh /path/to/int/test/root /path/to/env-file
#  ./paralleltestCreate.sh /path/to/int/test/root /path/to/env-file quick=1
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
commonFuncs="${geosChemDir}/test/shared/commonFunctionsForTests.sh"
. "${commonFuncs}"

# Echo header
printf "${SEP_MAJOR}\n"
printf "Creating GEOS-Chem Classic Parallelization Tests\n\n"
printf "GCClassic #${head_gcc}\n"
printf "GEOS_Chem #${head_gc}\n"
printf "HEMCO     #${head_hco}\n"
printf "${SEP_MAJOR}\n"

#=============================================================================
# Create integration test folder and subdirectories
#=============================================================================

# Create parallelization test root folder if it doesn't exist
ptRoot=$(absolute_path "${ptRoot}")
[[ ! -d "${ptRoot}" ]] && mkdir -p "${ptRoot}"

# Create local convenience variables
binDir="${ptRoot}/${BIN_DIR}"
buildDir="${ptRoot}/${BUILD_DIR}"
envDir="${ptRoot}/${ENV_DIR}"
execDir="${ptRoot}/${EXEC_DIR}"
logsDir="${ptRoot}/${LOGS_DIR}"
scriptsDir="${ptRoot}/${SCRIPTS_DIR}"
rundirsDir="${ptRoot}/${RUNDIRS_DIR}"

# Remove everything in the parallelization test root folder
cleanup_files "${ptRoot}"

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
printf "Linking to superproject      ${ptRoot}/CodeDir\n"
ln -s "${superProjectDir}" ${ptRoot}/CodeDir

#=============================================================================
# Copy files to the proper folders
#=============================================================================

printf "\nCopying run scripts to: ${ptRoot}\n"
cp -f ${envFile}                     ${envDir}/gcclassic.env
cp -f ${thisDir}/parallelTest*.sh    ${scriptsDir}
cp -f ${commonFuncs}                 ${scriptsDir}
cp -f ${thisDir}/README.md           ${scriptsDir}
cp -f ${thisDir}/README.testroot.md  ${ptRoot}/README.md

# Log file with echoback from rundir creation
log="${logsDir}/createParallelTests.log"

# Switch to folder where rundir creation scripts live
cd "${geosChemDir}/run/GCClassic"

#=============================================================================
# Create individual run directories: 4x5 - MERRA2 - 72L
#=============================================================================
printf "\nCreating new run directories:\n"

# 4x5 merra2 aerosol
create_rundir "2\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

#4x5 merra2 carbon
create_rundir "12\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

# 4x5 merra2 CH4
create_rundir "3\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

# 4x5 merra2 fullchem
create_rundir "1\n1\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

# DEBUG: Exit after creating a couple of rundirs if $quick is "yes"
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

# 4x5 merra2 metals"
create_rundir "11\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

#=============================================================================
# Create individual run directories: 4x5 and 47L (MERRA2)
#=============================================================================

# 4x5 merra2 fullchem_47L"
create_rundir "1\n1\n1\n1\n2\n${rundirsDir}\n\nn\n" "${log}"

#=============================================================================
# Nested-grid simulations
#=============================================================================

# 05x0625 merra2 CH4_47L_na"
create_rundir "3\n1\n3\n4\n2\n${rundirsDir}\n\nn\n" "${log}"

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
unset log
unset logsDir
unset ptRoot
unset rundirsDir
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
