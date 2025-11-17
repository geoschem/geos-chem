#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: redoIntegrationTestCompile.sh
#
# !DESCRIPTION: Manually resubmits an integration test compilation job.
#  Useful in case the original job died due to cluster issues, etc.
#\\
#\\
# !CALLING SEQUENCE:
#  cd /path/to/int/test/root/utils
#  ./redoIntegrationTestCompile.sh
#EOP
#------------------------------------------------------------------------------
#BOC

# Current directory
thisDir=$(realpath .)

# Load common functions
. "${thisDir}/../../shared/commonFunctionsForTests.sh"

# Do not let the tests proceed if a conda environment is activated,
# which could result in GEOS-Chem being linked against incorrect libraries!
activeCondaEnv=$(is_a_conda_env_activated)
if [[ "x${activeCondaEnv}" == "xtrue" ]]; then
    echo "ERROR: Tests cannot submitted when a conda environment is active!"
    echo "This may result in GEOS-Chem Classic being linked against the wrong libraries."
    echo "Please deactivate the environment and submit the tests again."
    exit 1
fi

# Integration test root dir is one dir higher
itRoot=$(realpath "${thisDir}/..")

# Change to the logs folder
cd "$itRoot/logs"

# Script to execute
script="${itRoot}/scripts/integrationTestCompile.sh"

# Make sure the script is found
if [[ ! -f "${script}" ]]; then
    echo "${script} was not found!  Exiting..."
    exit 1
fi

# Resubmit the compilation job to the scheduler
sbatch "${script}"
