#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: redoIntegrationTestExecute.sh
#
# !DESCRIPTION: Manually resubmits an integration test execution job.
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

# Integration test root dir is one dir higher
itRoot=$(realpath "${thisDir}/..")

# Change to the logs folder
cd "$itRoot/logs"

# Script to execute
script="${itRoot}/scripts/integrationTestExecute.sh"

# Make sure the script is found
if [[ ! -f "${script}" ]]; then
    echo "${script} was not found!  Exiting..."
    exit 1
fi

# Resubmit the compilation job to the scheduler
sbatch "${script}"
