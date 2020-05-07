#!/bin/bash

rm -f trac_avg.*
rm -f tracerinfo.dat
rm -f diaginfo.dat
rm -f cap_restart
rm -f gcchem*
rm -f *.rcx
rm -f *~
rm -f gchp.log
rm -f HEMCO.log
rm -f PET*.log
rm -f multirun.log
rm -f logfile.000000.out
rm -f slurm-*
rm -f 1
rm -f EGRESS

# Clean data too. Prompt user to confirm they want to do this.
# perhaps asking if they want to archive before deletion.
rm -f ./OutputDir/*.nc4
