#!/bin/bash

# setEnvironment: set local environment symbolic link gchp.env
#    
# Arguments: pass full path (no symlinks) to your customized GCHP env file
#
# Usage: ./setEnvironment /path/to/your/env/file
#
# E. Lundgren, 10/12/2018

# Make sure user passes environment file path
if [[ $# != 1 ]] ; then
 echo "Usage: ./setEnvironment /path/to/your/env/file"
 exit 1
fi

# Set symlink
if [[ -L gchp.env ]]; then
    unlink gchp.env
fi
if [[ -f "$1" ]]; then
  ln -s $1 gchp.env
  file gchp.env
else
  echo "Error: gchp.env target does not exist"
  exit 1
fi


exit 0
