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
#  ./intTestCreate.sh /path/to/integration/test/root/folder    # All tests
#  ./intTestCreate.sh /path/to/integration/test/root/folder 1  # debug, 2 tests
#
# !REMARKS:
#  Right now we pass values to the existing ./createRunDir.sh,
#  but will implement a more elegant solution later.#
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
ROOT=${1}
if [[ "x${ROOT}" == "x" ]]; then
    echo "ERROR: The root-level directory for tests has not been specified!"
    exit 1
fi

# Environment file
ENV_FILE=${2}
if [[ "x${ENV_FILE}" == "x" ]]; then
    echo "ERROR: The enviroment file (w/ module loads) has not been specified!"
    exit 1
fi
if [[ ! -f ${ENV_FILE} ]]; then
    echo "ERROR: The enviroment file is not a valid file!"
    exit 1
fi

# Debug switch
DEBUG=${3}
    
#=============================================================================
# Global variable and function definitions
#=============================================================================

# Current directory
TEST_DIR=$(pwd -P)
cd ${TEST_DIR}

# Top-level GEOS-Chem directory
cd ../..
GEOS_CHEM_DIR=$(pwd -P)

# GCClassic superproject directory
cd ../../
SUPER_PROJECT_DIR=$(pwd -P)
cd ${SUPER_PROJECT_DIR}

# Directory where the run creation scripts are found
RUN_DIR=${GEOS_CHEM_DIR}/run/GCClassic

# Load file with utility functions to setup configuration files
. ${GEOS_CHEM_DIR}/test/shared/commonFunctionsForTests.sh

# Get the absolute path of the root folder
ROOT=$(absolute_path ${ROOT})

# Log file
LOG=${ROOT}/logs/createIntTests.log

# Echo header
printf "${SEP_MAJOR}\n"
printf "Creating GEOS-Chem Classic Integration Tests\n"
printf "${SEP_MAJOR}\n"

#=============================================================================
# Initial setup of integration test directory
#=============================================================================

# Create integration test root folder if it doesn't exist
if [[ ! -d ${ROOT} ]]; then
    mkdir -p ${ROOT}
fi

# Remove run directories in the test folder
cleanup_files ${ROOT}

# Copying the run scripts to the Integration Test root folder
printf "\nCopying run scripts to: ${ROOT}\n"
cp ${ENV_FILE} ${ROOT}/gcclassic_env.sh
cp ${TEST_DIR}/intTest*.sh ${ROOT}
cp ${TEST_DIR}/commonFunctionsForTests.sh ${ROOT}

# Create a symbolic link to the code from the Integration Test root folder
ln -s ${SUPER_PROJECT_DIR} ${ROOT}/CodeDir

# Create log directory
if [[ !(-d ${ROOT}/logs) ]]; then
    printf "\nCreating log directory: ${ROOT}/logs\n"
    mkdir ${ROOT}/logs
fi

# Change to the directory where we will create the rundirs
printf "\nCreating new run directories:\n"

# Switch to folder where rundir creation scripts live
cd ${RUN_DIR}

#=============================================================================
# Create individual run directories: 2x25 - MERRA2 - 72L
#=============================================================================

DIR="merra2_2x25_CH4"
create_rundir "3\n1\n2\n1\n${ROOT}\n${DIR}\nn\n"          ${ROOT} ${DIR} ${LOG}

DIR="merra2_2x25_CO2"
create_rundir "4\n1\n2\n1\n${ROOT}\n${DIR}\nn\n"          ${ROOT} ${DIR} ${LOG}

# DEBUG: Exit after creating a couple of rundirs
# if the 2nd argument is passed and not a null string
if [[ "x${DEBUG}" != "x" ]]; then
    cd ${TEST_DIR}
    exit 0
fi

DIR="merra2_2x25_fullchem"
create_rundir "1\n1\n1\n1\n2\n1\n${ROOT}\n${DIR}\nn\n"    ${ROOT} ${DIR} ${LOG}

DIR="merra2_2x25_fullchem+aciduptake"
create_rundir "1\n1\n5\n1\n2\n1\n${ROOT}\n${DIR}\nn\n"    ${ROOT} ${DIR} ${LOG}

DIR="merra2_2x25_fullchem+APM"
create_rundir "1\n1\n7\n1\n2\n1\n${ROOT}\n${DIR}\nn\n"    ${ROOT} ${DIR} ${LOG}

DIR="merra2_2x25_fullchem+benchmark"
create_rundir "1\n1\n2\n1\n2\n1\n${ROOT}\n${DIR}\nn\n"    ${ROOT} ${DIR} ${LOG}

DIR="merra2_2x25_fullchem+complexSOA"
create_rundir "1\n1\n3\n1\n2\n1\n1\n${ROOT}\n${DIR}\nn\n" ${ROOT} ${DIR} ${LOG}

DIR="merra2_2x25_fullchem+complexSOA+SVPOA"
create_rundir "1\n1\n3\n2\n2\n1\n1\n${ROOT}\n${DIR}\nn\n" ${ROOT} ${DIR} ${LOG}

DIR="merra2_2x25_fullchem+marinePOA"
create_rundir "1\n1\n4\n1\n2\n1\n${ROOT}\n${DIR}\nn\n"    ${ROOT} ${DIR} ${LOG}

DIR="merra2_2x25_Hg"
create_rundir "5\n1\n2\n1\n${ROOT}\n${DIR}\nn\n"          ${ROOT} ${DIR} ${LOG}

DIR="merra2_2x25_POPs_BaP"
create_rundir "6\n1\n1\n2\n1\n${ROOT}\n${DIR}\nn\n"       ${ROOT} ${DIR} ${LOG}

DIR="merra2_2x25_tagCH4"
create_rundir "7\n1\n2\n1\n${ROOT}\n${DIR}\nn\n"          ${ROOT} ${DIR} ${LOG}

DIR="merra2_2x25_tagCO"
create_rundir "8\n1\n2\n1\n${ROOT}\n${DIR}\nn\n"          ${ROOT} ${DIR} ${LOG}

DIR="merra2_2x25_tagO3"
create_rundir "9\n1\n2\n1\n${ROOT}\n${DIR}\nn\n"          ${ROOT} ${DIR} ${LOG}

DIR="merra2_2x25_TransportTracers"
create_rundir "10\n1\n2\n1\n${ROOT}\n${DIR}\nn\n"         ${ROOT} ${DIR} ${LOG}

#=============================================================================
# Create individual run directories: 2x25 - GEOSFP - 72L
#=============================================================================

DIR="geosfp_2x25_CH4"
create_rundir "3\n2\n2\n1\n${ROOT}\n${DIR}\nn\n"          ${ROOT} ${DIR} ${LOG}

DIR="geosfp_2x25_CO2"
create_rundir "4\n2\n2\n1\n${ROOT}\n${DIR}\nn\n"          ${ROOT} ${DIR} ${LOG}

DIR="geosfp_2x25_fullchem"
create_rundir "1\n1\n1\n2\n2\n1\n${ROOT}\n${DIR}\nn\n"    ${ROOT} ${DIR} ${LOG}

DIR="geosfp_2x25_fullchem+aciduptake"
create_rundir "1\n1\n5\n2\n2\n1\n${ROOT}\n${DIR}\nn\n"    ${ROOT} ${DIR} ${LOG}

DIR="geosfp_2x25_fullchem+APM"
create_rundir "1\n1\n7\n2\n2\n1\n${ROOT}\n${DIR}\nn\n"    ${ROOT} ${DIR} ${LOG}

DIR="geosfp_2x25_fullchem+benchmark"
create_rundir "1\n1\n2\n2\n2\n1\n${ROOT}\n${DIR}\nn\n"    ${ROOT} ${DIR} ${LOG}

DIR="geosfp_2x25_fullchem+complexSOA"
create_rundir "1\n1\n3\n2\n2\n1\n1\n${ROOT}\n${DIR}\nn\n" ${ROOT} ${DIR} ${LOG}

DIR="geosfp_2x25_fullchem+complexSOA+SVPOA"
create_rundir "1\n1\n3\n2\n2\n1\n1\n${ROOT}\n${DIR}\nn\n" ${ROOT} ${DIR} ${LOG}

DIR="geosfp_2x25_fullchem+marinePOA"
create_rundir "1\n1\n4\n2\n2\n1\n${ROOT}\n${DIR}\nn\n"    ${ROOT} ${DIR} ${LOG}

DIR="geosfp_2x25_Hg"
create_rundir "5\n2\n2\n1\n${ROOT}\n${DIR}\nn\n"          ${ROOT} ${DIR} ${LOG}

DIR="geosfp_2x25_POPs_BaP"
create_rundir "6\n1\n2\n2\n1\n${ROOT}\n${DIR}\nn\n"       ${ROOT} ${DIR} ${LOG}

DIR="geosfp_2x25_tagCH4"
create_rundir "7\n2\n2\n1\n${ROOT}\n${DIR}\nn\n"          ${ROOT} ${DIR} ${LOG}

DIR="geosfp_2x25_tagCO"
create_rundir "8\n2\n2\n1\n${ROOT}\n${DIR}\nn\n"          ${ROOT} ${DIR} ${LOG}

DIR="geosfp_2x25_tagO3"
create_rundir "9\n2\n2\n1\n${ROOT}\n${DIR}\nn\n"          ${ROOT} ${DIR} ${LOG}

DIR="geosfp_2x25_TransportTracers"
create_rundir "10\n2\n2\n1\n${ROOT}\n${DIR}\nn\n"         ${ROOT} ${DIR} ${LOG}

#=============================================================================
# Create individual run directories: 4x5 - MERRA2 - 72L
#=============================================================================

DIR="merra2_4x5_CH4"
create_rundir "3\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"          ${ROOT} ${DIR} ${LOG}

DIR="merra2_4x5_fullchem"
create_rundir "1\n1\n1\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"    ${ROOT} ${DIR} ${LOG}

DIR="merra2_4x5_fullchem+aciduptake"
create_rundir "1\n1\n5\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"    ${ROOT} ${DIR} ${LOG}

DIR="merra2_4x5_fullchem+APM"
create_rundir "1\n1\n7\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"    ${ROOT} ${DIR} ${LOG}

DIR="merra2_4x5_fullchem+benchmark"
create_rundir "1\n1\n2\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"    ${ROOT} ${DIR} ${LOG}

DIR="merra2_4x5_fullchem+complexSOA"
create_rundir "1\n1\n3\n1\n1\n1\n1\n${ROOT}\n${DIR}\nn\n" ${ROOT} ${DIR} ${LOG}

DIR="merra2_4x5_fullchem+complexSOA+SVPOA"
create_rundir "1\n1\n3\n2\n1\n1\n1\n${ROOT}\n${DIR}\nn\n" ${ROOT} ${DIR} ${LOG}

DIR="merra2_4x5_fullchem+marinePOA"
create_rundir "1\n1\n4\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"    ${ROOT} ${DIR} ${LOG}

DIR="merra2_4x5_fullchem+RRTMG"
create_rundir "1\n1\n8\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"    ${ROOT} ${DIR} ${LOG}

DIR="merra2_4x5_fullchem+TOMAS15"
create_rundir "1\n1\n6\n1\n1\n1\n1\n${ROOT}\n${DIR}\nn\n" ${ROOT} ${DIR} ${LOG}

DIR="merra2_4x5_Hg"
create_rundir "5\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"          ${ROOT} ${DIR} ${LOG}

DIR="merra2_4x5_POPs_BaP"
create_rundir "6\n1\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"       ${ROOT} ${DIR} ${LOG}

DIR="merra2_4x5_tagCH4"
create_rundir "7\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"          ${ROOT} ${DIR} ${LOG}

DIR="merra2_4x5_tagCO"
create_rundir "8\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"          ${ROOT} ${DIR} ${LOG}

DIR="merra2_4x5_tagO3"
create_rundir "9\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"          ${ROOT} ${DIR} ${LOG}

DIR="merra2_4x5_TransportTracers"
create_rundir "10\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"         ${ROOT} ${DIR} ${LOG}

#=============================================================================
# Create individual run directories: 4x5 - GEOSFP - 72L
#=============================================================================

DIR="geosfp_4x5_CH4"
create_rundir "3\n2\n1\n1\n${ROOT}\n${DIR}\nn\n"          ${ROOT} ${DIR} ${LOG}

DIR="geosfp_4x5_fullchem"
create_rundir "1\n1\n1\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"    ${ROOT} ${DIR} ${LOG}

DIR="geosfp_4x5_fullchem+aciduptake"
create_rundir "1\n1\n5\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"    ${ROOT} ${DIR} ${LOG}

DIR="geosfp_4x5_fullchem+APM"
create_rundir "1\n1\n7\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"    ${ROOT} ${DIR} ${LOG}

DIR="geosfp_4x5_fullchem+benchmark"
create_rundir "1\n1\n2\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"    ${ROOT} ${DIR} ${LOG}

DIR="geosfp_4x5_fullchem+complexSOA"
create_rundir "1\n1\n3\n1\n1\n1\n1\n${ROOT}\n${DIR}\nn\n" ${ROOT} ${DIR} ${LOG}

DIR="geosfp_4x5_fullchem+complexSOA+SVPOA"
create_rundir "1\n1\n3\n2\n1\n1\n1\n${ROOT}\n${DIR}\nn\n" ${ROOT} ${DIR} ${LOG}

DIR="geosfp_4x5_fullchem+marinePOA"
create_rundir "1\n1\n4\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"    ${ROOT} ${DIR} ${LOG}

DIR="geosfp_4x5_fullchem+RRTMG"
create_rundir "1\n1\n8\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"    ${ROOT} ${DIR} ${LOG}

DIR="geosfp_4x5_fullchem+TOMAS15"
create_rundir "1\n1\n6\n1\n1\n1\n1\n${ROOT}\n${DIR}\nn\n" ${ROOT} ${DIR} ${LOG}

DIR="geosfp_4x5_Hg"
create_rundir "5\n2\n1\n1\n${ROOT}\n${DIR}\nn\n"          ${ROOT} ${DIR} ${LOG}

DIR="geosfp_4x5_POPs_BaP"
create_rundir "6\n2\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"       ${ROOT} ${DIR} ${LOG}

DIR="geosfp_4x5_tagCH4"
create_rundir "7\n2\n1\n1\n${ROOT}\n${DIR}\nn\n"          ${ROOT} ${DIR} ${LOG}

DIR="geosfp_4x5_tagCO"
create_rundir "8\n2\n1\n1\n${ROOT}\n${DIR}\nn\n"          ${ROOT} ${DIR} ${LOG}

DIR="geosfp_4x5_tagO3"
create_rundir "9\n2\n1\n1\n${ROOT}\n${DIR}\nn\n"          ${ROOT} ${DIR} ${LOG}

DIR="geosfp_4x5_TransportTracers"
create_rundir "10\n2\n1\n1\n${ROOT}\n${DIR}\nn\n"         ${ROOT} ${DIR} ${LOG}

# Comment these out for now (bmy, 11/12/20)
##=============================================================================
## Nested-grid simulations 
##=============================================================================
#
#DIR="merra2_05x0625_CH4_na_47L"
#create_rundir "3\n1\n3\n4\n2\n${ROOT}\n${DIR}\nn\n"       ${ROOT} ${DIR} ${LOG}
#
#DIR="geosfp_025x03125_CH4_na_47L"
#create_rundir "3\n2\n4\n4\n2\n${ROOT}\n${DIR}\nn\n"       ${ROOT} ${DIR} ${LOG}
#
#DIR="merra2_05x0625_fullchem_as_47L"
#create_rundir "1\n2\n1\n1\n3\n2\n2\n${ROOT}\n${DIR}\nn\n" ${ROOT} ${DIR} ${LOG}
#
#DIR="geosfp_025x03125_fullchem_na_47L"
#create_rundir "1\n2\n1\n2\n4\n4\n2\n${ROOT}\n${DIR}\nn\n" ${ROOT} ${DIR} ${LOG}

#=============================================================================
# Cleanup and quit
#=============================================================================

# Switch back to the present directory
cd ${TEST_DIR}

# Free local variables
unset ROOT
unset TEST_DIR
unset GEOS_CHEM_DIR
unset SUPER_PROJECT_DIR
unset RUN_DIR
unset LOG
unset DIR

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
