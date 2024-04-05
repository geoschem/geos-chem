#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: redoParallelTestCompile.sh
#
# !DESCRIPTION: Manually resubmits an parallel test compilation job.
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
ptRoot=$(realpath "${thisDir}/..")

# Change to the logs folder
cd "$itRoot/logs"

# Script to execute
script="${ptRoot}/scripts/parallelTestCompile.sh"

# Make sure the script is found
if [[ ! -f "${script}" ]]; then
    echo "${script} was not found!  Exiting..."
    exit 1
fi

# Resubmit the compilation job to the scheduler
sbatch "${script}"
