#!/bin/bash

# Script to archive files after a run. 
# 
# Argument: archive directory name (can be non-existent)
#
# Example usage: ./archiveRun.sh 1mon_24hrDiag
#
# The output data (OutputDir/*.nc4) is moved but everything else is copied, 
# including log files (*.log, slurm-*), config files (*.rc, input.geos), 
# run files (*.run), and restarts (GEOSChem.Restart.*, HEMCO_restart.*). 
# Files are stored in subdirectories within the archive directory.
#
# NOTE: Clean the run directory AFTER archiving with './cleanupRunDir.sh'
# if you plan on doing another run. Otherwise previous run files will also
# be archived if this script is called again.

# Initial version: Lizzie Lundgren - 7/12/2018

# Customize this script as needed to best fit your workflow.

# Check that directory name passed
if [[ $# == 1 ]]; then
    archivedir=$1
else
   echo "Usage: ./archiveRun.sh {ArchiveDirName}"
   exit 
fi

# Check that directory does not already exist
if [ -d "${archivedir}" ]; then
   echo "Warning: Directory ${archivedir} already exists."
   echo "Remove or rename that directory, or choose a different name."
   exit 1
fi

# Function to move files and subdirs in directory except if string match
# ( arg1 : source, arg2 : target, arg3 : exclude string )
movefiles () {
   numMoved=0
   for item in $1/*; do
      if [[ $(basename $item) == $3 ]]; then
         continue
      elif [[ -e $item ]]; then
         if [[ -d $item ]]; then
            echo "   -> $2/$(basename $item)/"
         else
            echo "   -> $2/$(basename $item)"
         fi
         mv $item $2
         numMoved=$numMoved+1
      fi
   done
   if [[ $numMoved == "0" ]]; then
      echo "   Warning: No files to move from $1" 
   fi
}

# Function to copy all files matching string (arg2) to directory (arg1)
# ( arg1 : source, arg2 : target )
copyfiles () {
   for file in $1; do
      if [ -e $file ]; then
         echo "   -> $2/$file"
         cp -t $2 $file
      else
         echo "   Warning: $file not found"
      fi
   done
}

# Make Archive directory
echo "Archiving files to directory $1"
mkdir -p ${archivedir}
mkdir -p ${archivedir}/OutputDir
mkdir -p ${archivedir}/BenchmarkResults
mkdir -p ${archivedir}/Logs
mkdir -p ${archivedir}/Config
mkdir -p ${archivedir}/Restarts
mkdir -p ${archivedir}/Build

# Move large files rather than copy (except initial restart)
echo "Moving files and directories..."
movefiles "BenchmarkResults" ${archivedir}/BenchmarkResults
movefiles "OutputDir" ${archivedir}/OutputDir FILLER

# Copy everything else
echo "Copying files..."
copyfiles "input.geos"         ${archivedir}/Config
copyfiles "rundir.version"     ${archivedir}/Config
copyfiles "*.rc"               ${archivedir}/Config
copyfiles "*.run"              ${archivedir}/Config
copyfiles "*.env"              ${archivedir}/Config
copyfiles "*.log"              ${archivedir}/Logs
copyfiles "slurm-*"            ${archivedir}/Logs
copyfiles "GEOSChem.Restart.*" ${archivedir}/Restarts
copyfiles "HEMCO_restart.*"    ${archivedir}/Restarts
copyfiles "build_info/*"       ${archivedir}/Build

printf "Complete!\n"

exit 0
