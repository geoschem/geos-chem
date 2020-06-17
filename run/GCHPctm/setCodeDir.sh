#!/bin/bash

# setCodeDir: bash script to set local CodeDir symbolic link for GEOS-Chem source code
#    
# Arguments: pass full path (no symlinks) to GEOS-Chem source code directory
#
# Usage: ./setCodeDir.sh /path/to/your/code
#
# E. Lundgren, 8/17/2017

# Prompt the user for source code directory and set symlink
if [[ -L CodeDir ]]; then
    unlink CodeDir
fi
if [[ -d "$1" ]]; then
  ln -s $1 CodeDir
  file CodeDir
  printf "\nPrinting source code repository information...\n\n"
  make printbuildinfo
  printf "\nType 'make printbuildinfo' at any time in this run directory to view source code status.\n"
else
  echo "Error: Source code directory does not exist"
  exit 1
fi


exit 0
