#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: createIntTests.sh
#
# !DESCRIPTION: Creates integration test run directories in a user-specified
#  root folder, and copies a run script there.
#\\
#\\
# !CALLING SEQUENCE:
#  ./createIntTests /path/to/integration/test/root/folder
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
# Check arguments
#=============================================================================
if [[ "x${1}" == "x" ]]; then
    echo "ERROR: The root-level directory for tests has not been specified!"
    exit 1
fi

if [[ "x${2}" == "x" ]]; then
    echo "ERROR: The environment file has not been specified!"
    exit 1
fi

#=============================================================================
# Global variable and function definitions
#=============================================================================

# Integration test root folder (create if it doesn't exist)
ROOT=${1}
if [[ !(-d ${ROOT}) ]]; then
    mkdir -p ${ROOT}
fi

# Environment file
ENV_FILE=${2}

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
RUN_DIR=${GEOS_CHEM_DIR}/run/GCHPctm

# Load file with utility functions to setup configuration files
. ${GEOS_CHEM_DIR}/test/shared/commonFunctionsForTests.sh

# Get the absolute path of the root folder
ROOT=$(absolute_path ${ROOT})

# Log file
LOG=${ROOT}/logs/createIntTests.log

# Echo header
printf "${SEP_MAJOR}\n"
printf "Creating GCHPctm Integration Tests\n"
printf "${SEP_MAJOR}\n"

#=============================================================================
# Initial setup of integration test directory
#=============================================================================

# Remove run directories in the test folder
cleanup_files ${ROOT}

# Copying the run scripts to the root test folder
printf "\nCopying run scripts to: ${ROOT}\n"
cp ${TEST_DIR}/intTest*.sh ${ROOT}
cp ${TEST_DIR}/commonFunctionsForTests.sh ${ROOT}

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
# Create individual run directories
#=============================================================================

DIR="merra2_TransportTracers"
create_rundir "2\n1\n${ROOT}\n${DIR}\nn\n" ${ROOT} ${DIR} ${LOG}
cp ${ENV_FILE} ${ROOT}/${DIR}/gchp.env
cp ${TEST_DIR}/gchp.slurm.sh ${ROOT}/${DIR}/gchp.slurm.sh
cp ${TEST_DIR}/gchp.lsf.sh   ${ROOT}/${DIR}/gchp.lsf.sh

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
