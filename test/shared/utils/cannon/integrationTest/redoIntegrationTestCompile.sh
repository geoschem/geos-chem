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

# Throw error if there is a conda environment active with netCDF,
# which can cause the code to be linked against the wrong netCDF version.
if [[ -n "$CONDA_DEFAULT_ENV" || \
      "$(which nc-config 2>/dev/null)" == *conda* ]]; then
   echo "ERROR: Conda netCDF detected. Run 'conda deactivate' first."
   exit 1
fi

# Current directory
thisDir=$(realpath .)

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
