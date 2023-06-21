#!/bin/bash

rm -fv trac_avg.*
rm -fv tracerinfo.dat
rm -fv diaginfo.dat
rm -fv gcchem*
rm -fv *.rcx
rm -fv *~
rm -fv gchp.log
rm -fv gchp.*.log
rm -fv HEMCO.log
rm -fv PET*.log
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

# Clean data too.  If an argument is passed, then prompt user to confirm
# perhaps asking if they want to archive before deletion.
if [[ "x${1}" == "x" ]]; then
    rm -Iv ./OutputDir/*.nc*     # Get confirmation from user
else
    rm -fv ./OutputDir/*.nc*     # Skip confirmation from user
fi

# Give instruction to reset start date if using GCHP
echo "Reset simulation start date in cap_restart if using GCHP"
