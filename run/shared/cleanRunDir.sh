#!/bin/bash

#============================================================================
# cleanRunDir.sh: Removes files created by GEOS-Chem from a run directory
#
# Usage:
# ------
# $ ./cleanRunDir.sh          # Removes model output files in the run
#                             # directory.  Also prompts the user before
#                             # removing diagnostic output files from
#                             # from OutputDir/.
#
# $ ./cleanRunDir.sh --force  # Removes model output files in the run
#                             # directory, but will remove diagnostic
#                             # output files without prompting first.
#                             # USE WITH CAUTION!
#============================================================================

# Clean model output files in the run directory
rm -fv gcchem*
rm -fv *.rcx
rm -fv *~
rm -fv gchp.log
rm -fv gchp.*.log
rm -fv HEMCO.log
rm -fv PET*.log
rm -fv ESMF_LogFile
rm -fv multirun.log
rm -fv warnings_and_errors.log
rm -fv GC*.log
rm -fv log.dryrun*
rm -fv logfile.000000.out
rm -fv slurm-*
rm -fv 1
rm -fv EGRESS
rm -fv core.*
rm -fv PET*.ESMF_LogFile
rm -fv allPEs.log
rm -fv fort.*

#----------------------------------------------------------------------------
# Clean data files in OutputDir.
# These are netCDF files (*.nc) and KPP standalone interface files (*.txt).
#----------------------------------------------------------------------------
if [[ "x${1}" == "x" ]]; then      # User confirmation required
    rm -Iv ./OutputDir/*.nc*
    rm -Iv ./OutputDir/*.txt
else                               # User Confirmation not required
    rm -fv ./OutputDir/*.nc*
    rm -fv ./OutputDir/*.txt*
fi

#---------------------------------------------------------------------------
# Give instruction to reset start date if using GCHP
#---------------------------------------------------------------------------
echo "Reset simulation start date in cap_restart if using GCHP"
