#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: intTests_lsf.sh
#
# !DESCRIPTION: Runs integration tests on the various GEOS-Chem Classic
#  run directories (using the LSF scheduler).  Compilation tests and
#  execution tests are submitted as separate SLURM jobs, using dependencies.
#\\
#\\
# !CALLING SEQUENCE:
#  ./intTests_lsf.sh /integration/test/root/folder ENV-FILE          or
#  ./intTests_lsf.sh /integration/test/root/folder ENV-FILE short=1
#
# !REVISION HISTORY:
#  03 Nov 2020 - R. Yantosca - Initial version
#  See the subsequent Git history with the gitk browser!
#EOP
#------------------------------------------------------------------------------
#BOC

#=============================================================================
# Check arguments and define variables
#=============================================================================

# 1st argument: Root directory for integration tests
INT_TEST_ROOT=${1}
if [[ "x${INT_TEST_ROOT}" == "x" ]]; then
    echo "ERROR: The root-level directory for tests has not been specified!"
    exit
fi

# 2nd argument: Environment file
ENV_FILE=${2}
if [[ "x${ENV_FILE}" == "x" ]]; then
    echo "ERROR: The enviroment file (w/ module loads) has not been specified!"
    exit
fi

# 3rd argument: Run a short integration test (for development)?
SHORT=${3}

#=============================================================================
# Load file with utility functions to setup configuration files
#=============================================================================

# Current directory
THIS_DIR=$(pwd -P)

# Load common functions
. ${THIS_DIR}/commonFunctionsForTests.sh

#=============================================================================
# Create integration test directories in the root folder
#=============================================================================

# Convert integration test root folder to an absolute path
INT_TEST_ROOT=$(absolute_path ${INT_TEST_ROOT})

# Create GEOS-Chem run directories in the integration test root folder
./intTestCreate.sh ${INT_TEST_ROOT} ${ENV_FILE} ${SHORT}
if [[ $? -ne 0 ]]; then
   exit 0
fi

# Change to the integration test root folder
if [[ -d ${INT_TEST_ROOT} ]]; then
    cd ${INT_TEST_ROOT}
else
    echo "ERROR: ${INT_TEST_ROOT} is not a valid directory!  Exiting..."
    exit 1
fi

#=============================================================================
# Submit compilation tests script
#
#### Liam: Modify code to return job id

#=============================================================================
OUTPUT=$(bsub intTestCompile_slurm.sh)
OUTPUT=($OUTPUT)
CMP_ID=${OUTPUT[3]}

#=============================================================================
# Submit execution tests script as a job dependency
#
#### Liam: Modify code to submit LSF job depending on prior job
#=============================================================================
OUTPUT=$(bsub --dependency... intTestExecute_slurm.sh)
OUTPUT=($OUTPUT)
EXE_ID=${OUTPUT[3]}

# Change back to this directory
cd ${THIS_DIR}

#=============================================================================
# Cleanup and quit
#=============================================================================
echo ""
echo "Compilation tests submitted as LSF job ${CMP_ID}"
echo "Execution   tests submitted as LSF job ${EXE_ID}"

# Free local variables
unset ENV_FILE
unset INT_TEST_ROOT
unset JOBID
unset OUTPUT
unset SHORT
unset THIS_DIR

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
