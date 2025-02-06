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
#  ./integrationTestCreate.sh /path/to/int/test/root /path/to/env-file [yes|no]
#EOP
#------------------------------------------------------------------------------
#BOC

#=============================================================================
# Load common functions
#=============================================================================

# Current directory
thisDir=$(pwd -P)
cd "${thisDir}"

# Path to the test/shared folder in source code
sharedDir=$(realpath "${thisDir}/../../shared")

# Source the script containing utility functions and variables
commonFuncs="${sharedDir}/commonFunctionsForTests.sh"
. "${commonFuncs}"

#=============================================================================
# Parse input arguments
#=============================================================================

# Integration test root folder
itRoot="${1}"
if [[ "X${itRoot}" == "X" ]]; then
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

# Run a compile-only integration test?
testsToRun="${3}"

# Run a short integration test?
quick="${4}"

#=============================================================================
# Global variable and function definitions
#=============================================================================

# GCClassic superproject directory (absolute paths)
cd ../../../../../
superProjectDir=$(pwd -P)
cd ${superProjectDir}

# GEOS-Chem and HEMCO submodule directories
geosChemDir="${superProjectDir}/src/GEOS-Chem"

# Echo header
printf "${SEP_MAJOR}\n"
printf "Creating GEOS-Chem Classic Integration Tests\n\n"
print_submodule_head_commits "10" "${superProjectDir}" ""
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
utilsDir="${itRoot}/${UTILS_DIR}"

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

# Subdir for env files (for Harvard Cannon only)
if [[ "X${site}" == "XCANNON" ]]; then
    printf "Creating env files directory ${envDir}\n"
    mkdir -p "${envDir}"
fi

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

printf "\nCopying run scripts to:      ${scriptsDir}\n"
cp -f ${thisDir}/integration*.sh     ${scriptsDir}
cp -f ${commonFuncs}                 ${scriptsDir}
cp -f ${thisDir}/README.md           ${scriptsDir}
cp -f ${thisDir}/README.testroot.md  ${itRoot}/README.md

if [[ "X${site}" == "XCANNON" ]]; then

    # Copy Cannon environment file
    cp -f  ${envFile} ${envDir}/gcclassic.env

    # Copy Cannon utility scripts
    printf "Copying utility scripts to   ${utilsDir}\n"
    cp -fR ${sharedDir}/utils/cannon/integrationTest  ${utilsDir}

elif [[ "X${site}" == "XCOMPUTE1" ]]; then

    # Copy Compute1 utility scripts
    printf "Copying utility scripts to   ${utilsDir}\n"
    cp -fR ${sharedDir}/utils/compute1/integrationTest  ${utilsDir}

    # Force scripts to be executable (Compute1 resets permissions)
    chmod 755 -R ${scriptsDir}
    chmod 755 -R ${utilsDir}

fi

# Log file with echoback from rundir creation
log="${logsDir}/createIntegrationTests.log"

#=============================================================================
# Don't create run directories for compile-only tests.
#=============================================================================
if [[ "X${testsToRun}" == "XALL" ]]; then

    # Switch to folder where rundir creation scripts live
    cd "${geosChemDir}/run/GCClassic"

    #=========================================================================
    # Create individual run directories: 4x5 - MERRA2 - 72L
    #=========================================================================
    printf "\nCreating new run directories:\n"

    # 4x5 merra2 aerosol
    create_rundir "2\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

    # 4x5 merra2 carbon
    create_rundir "3\n1\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

    # 4x5 merra2 carbon CH4 only
    create_rundir "3\n2\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

    # 4x5 merra2 carbon CO2 only
    create_rundir "3\n3\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

    # 4x5 merra2 carbon CO only
    create_rundir "3\n4\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

    # 4x5 merra2 carbon OCS only
    create_rundir "3\n5\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

    # 4x5 merra2 fullchem
    create_rundir "1\n1\n1\n1\n1\n${rundirsDir}\n\nn\nn\n" "${log}"

    # Exit after creating a couple of rundirsDirs if $quick is "yes"
    if [[ "X${quick}" == "XYES" ]]; then
        cd ${thisDir}
        exit 0
    fi

    # 4x5 merra2 fullchem_LuoWd
    dir="gc_4x5_merra2_fullchem_LuoWd"
    create_rundir "1\n1\n1\n1\n1\n${rundirsDir}\n${dir}\nn\nn\n" "${log}"

    # 4x5 merra2 fullchem_aciduptake
    create_rundir "1\n5\n1\n1\n1\n${rundirsDir}\n\nn\nn\n" "${log}"

    # 4x5 merra2 fullchem_APM
    create_rundir "1\n7\n1\n1\n1\n${rundirsDir}\n\nn\nn\n" "${log}"

    # 4x5 merra2 fullchem_benchmark
    create_rundir "1\n2\n1\n1\n1\n${rundirsDir}\n\nn\nn\n" "${log}"

    # 4x5 merra2 fullchem_complexSOA
    create_rundir "1\n3\n1\n1\n1\n1\n${rundirsDir}\n\nn\nn\n" "${log}"

    # 4x5 merra2 fullchem_complexSOA_SVPOA
    create_rundir "1\n3\n2\n1\n1\n1\n${rundirsDir}\n\nn\nn\n" "${log}"

    # 4x5 merra2 fullchem_marinePOA
    create_rundir "1\n4\n1\n1\n1\n${rundirsDir}\n\nn\nn\n" "${log}"

    # 4x5 merra2 fullchem_RRTMG
    create_rundir "1\n8\n1\n1\n1\n${rundirsDir}\n\nn\nn\n" "${log}"

    # 4x5 merra2 fullchem_TOMAS15_47L
    create_rundir "1\n6\n1\n1\n1\n2\n${rundirsDir}\n\nn\nn\n" "${log}"

    # 4x5 merra2 Hg
    create_rundir "4\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

    # 4x5 merra2 POPs_BaP
    create_rundir "5\n1\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

    # 4x5 merra2 tagO3
    create_rundir "6\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

    # 4x5 merra2 TransportTracers
    create_rundir "7\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

    # 4x5 merra2 TransportTracers_LuoWd
    dir="gc_4x5_merra2_TransportTracers_LuoWd"
    create_rundir "7\n1\n1\n1\n${rundirsDir}\n${dir}\nn\n" "${log}"

    # 4x5 merra2 metals
    create_rundir "8\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

    # 4x5 merra2 CH4
    create_rundir "9\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

    # 4x5 merra2 CO2
    create_rundir "10\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

    # 4x5 merra2 tagCO
    create_rundir "11\n1\n1\n1\n${rundirsDir}\n\nn\n" "${log}"

    #=========================================================================
    # Create individual run directories: 4x5 and 47L (MERRA2)
    #=========================================================================

    # 4x5 merra2 fullchem_47L
    create_rundir "1\n1\n1\n1\n2\n${rundirsDir}\n\nn\nn\n" "${log}"

    #=========================================================================
    # GCAP 2.0 simulations
    #=========================================================================

    # 2x2.5 ModelE2.1 fullchem (scenario SSP2-4.5, option 6)
    create_rundir "1\n1\n4\n6\n2\n1\n${rundirsDir}\n\nn\nn\n" "${log}"

    #=========================================================================
    # Nested-grid simulations
    #=========================================================================

    # 05x0625 merra2 CH4_47L_na
    create_rundir "9\n1\n3\n4\n2\n${rundirsDir}\n\nn\n" "${log}"

    # 05x0625 merra2 fullchem_47L_na
    create_rundir "1\n1\n1\n3\n4\n2\n${rundirsDir}\n\nn\nn\n" "${log}"

    #=========================================================================
    # Simulation with all diagnostics on
    #==========================================================================

    # Configuration files
    allDiagsDir="gc_4x5_merra2_fullchem_alldiags"
    extDataDir=$(grep "GC_DATA_ROOT" "${HOME}/.geoschem/config")
    extDataDir=${extDataDir/export GC_DATA_ROOT\=/}
    pfDat="${geosChemDir}/test/shared/alldiags/Planeflight.dat.20190701"
    obsPk="${extDataDir}/Data_for_Int_Tests/obspack_input_for_testing.20190701.nc"
    # Copy the fullchem_benchmark rundir to fullchem_alldiags
    echo "... ${itRoot}/rundirs/${allDiagsDir}"
    cd "${rundirsDir}"
    cp -r "gc_4x5_merra2_fullchem_benchmark" "${allDiagsDir}"
    cd "${allDiagsDir}"

    # Turn on all collections except RRTMG and Tomas collections (which
    # Make sure to activate these in the RRTMG and TOMAS integration tests.
    # Also note; there is a floating point error in the UVFlux diagnostic,
    # so temporarily comment that out.
    sed_ie "s|#'|'|"               "HISTORY.rc"
    sed_ie "s|'RRTMG'|#'RRTMG'|"   "HISTORY.rc"
    sed_ie "s|'Tomas'|#'Tomas'|"   "HISTORY.rc"
    sed_ie "s|'DynHeat|#'DynHeat|" "HISTORY.rc"
    sed_ie "s|'UVFlux'|#'UVFlux'|" "HISTORY.rc"

    # Activate the planeflight diagnostic
    cp -r "${pfDat}" .
    toggle_geoschem_config_option "geoschem_config.yml" "planeflight" "true "

    # Activate the ObsPack diagnostic
    cp -r "${obsPk}" .
    toggle_geoschem_config_option "geoschem_config.yml" "obspack"     "true "

    # Switch back to the present directory
    cd "${thisDir}"

fi

#=============================================================================
# Cleanup and quit
#=============================================================================

# Free local variables
unset allDiagsDir
unset binDir
unset buildDir
unset commonFuncs
unset dir
unset envDir
unset extDataDir
unset geosChemDir
unset itRoot
unset log
unset logsDir
unset pfDat
unset obsPk
unset rundirsDir
unset superProjectDir
unset scriptsDir
unset sharedDir
unset thisDir
unset utilsDir

# Free imported variables
unset FILL
unset SEP_MAJOR
unset SEP_MINOR
unset CMP_PASS_STR
unset CMP_FAIL_STR
unset EXE_PASS_STR
unset EXE_FAIL_STR
#EOC
