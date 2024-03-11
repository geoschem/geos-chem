#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: integrationTestCreate.sh
#
# !DESCRIPTION: Creates GCHP integration test run directories in a
#  user-specified root folder, and copies a run script there.
#\\
#\\
# !CALLING SEQUENCE:
#  ./integrationTestCreate.sh /path/to/root /path/to/env-file tests-to-run
#  ./integrationTestCreate.sh /path/to/root /path/to/env-file tests-to-run quick=1
#EOP
#------------------------------------------------------------------------------
#BOC

#=============================================================================
# Load common functions
#=============================================================================

# Current directory
thisDir=$(pwd -P)
cd "${thisDir}"

# Source the script containing utility functions and variables
. "${thisDir}../../shared/commonFunctionsForTests.sh"

#=============================================================================
# Parse input arguments
#=============================================================================

# Integration test root folder
itRoot="${1}"
if [[ "x${itRoot}" == "x" ]]; then
    echo "ERROR: The root-level directory for tests has not been specified!"
    exit 1
fi

# Environment file (for Harvard Cannon only)
site=$(get_site_name)
envFile="${2}"
if [[ "X${site}" == "XCANNON" ]]; then
    [[ "X${envFile}" == "X" ]] && envFile=$(get_default_gcc_env_file)
    if [[ ! -f ${envFile} ]]; then
	echo "ERROR: The enviroment file is not a valid file!"
	exit 1
    fi
fi

# Type of tests to run?
testsToRun="${3}"

# Run a short integration test?
quick="${4}"

#=============================================================================
# Global variable and function definitions
#=============================================================================

# GCClassic superproject directory
cd ../../../../../../../
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

# Echo header
printf "${SEP_MAJOR}\n"
printf "Creating GCHP Integration Tests\n\n"
printf "GCHP      #${head_gchp}\n"
printf "GEOS_Chem #${head_gc}\n"
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
if [[ "x${testsToRun}" == "xALL" ]]; then
    printf "Creating rundirs directory   ${rundirsDir}\n"
    mkdir -p "${rundirsDir}"
fi
    
# Create a symbolic link to the code from the Integration Test root folder
printf "Linking to superproject      ${itRoot}/CodeDir\n"
ln -s "${superProjectDir}" ${itRoot}/CodeDir

#=============================================================================
# Copy files to the proper folders
#=============================================================================

printf "\nCopying run scripts to: ${itRoot}/${SCRIPTS_DIR}\n"
cp -f ${envFile}                     ${envDir}/gchp.env
cp -f ${thisDir}/integration*.sh     ${scriptsDir}
cp -f ${commonFuncs}                 ${scriptsDir}
cp -f ${thisDir}/README.md           ${scriptsDir}
cp -f ${thisDir}/README.testroot.md  ${itRoot}/README.md

# This is necessary on Compute1 to make all scripts executable
chmod 755 -R ${scriptsDir}

# Log file with echoback from rundir creation
log="${logsDir}/createIntegrationTests.log"

#=============================================================================
# Don't create run directories for compile-only tests.
#=============================================================================
if [[ "X${testsToRun}" == "XALL" ]]; then

    # Switch to folder where rundir creation scripts live
    cd "${geosChemDir}/run/GCHP"

    #=========================================================================
    # Create the GCHP run directories
    #=========================================================================
    printf "\nCreating new run directories:\n"

    # c24 geosfp TransportTracers
    create_rundir "2\n1\n${rundirsDir}\n\nn\n" "${log}"

    # c24 merra2 fullchem tagO3
    create_rundir "4\n1\n${rundirsDir}\n\nn\n" "${log}"

    # Placeholder for carbon simulation
    # c24 merra2 carbon
    #create_rundir "12\n1\n${rundirsDir}\n\nn\n" "${log}"

    # DEBUG: Exit after creating a couple of rundirs if $quick is "yes"
    if [[ "X${quick}" == "XYES" ]]; then
        cd ${thisDir}
        exit 0
    fi

    # c24 merra2 fullchem_standard
    create_rundir "1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

    # c24 merra2 fullchem_benchmark
    create_rundir "1\n2\n1\n${rundirsDir}\n\nn\n" "${log}"

    # c24 merra2 fullchem_RRTMG
    create_rundir "1\n8\n1\n${rundirsDir}\n\nn\n" "${log}"

    # c24 merra2 fullchem_TOMAS15
    create_rundir "1\n6\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

    # Switch back to the present directory
    cd "${thisDir}"
fi
    
#=============================================================================
# Cleanup and quit
#=============================================================================

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
