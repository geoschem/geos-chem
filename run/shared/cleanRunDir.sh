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

# Clean data too. Prompt user to confirm they want to do this.
# perhaps asking if they want to archive before deletion.
rm -Iv ./OutputDir/*.nc*
