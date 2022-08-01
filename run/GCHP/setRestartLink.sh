#!/bin/bash

# setRestartLink.sh
#
# DESCRIPTION:
#    Sets local restart symbolic link gchp_restart.nc4 to point to
#    ./Restart/GEOSChem.Restart.YYYYMMDD_HHmmz.cN.nc4, where YYYYMMDD and
#    HHmm are the simulation start date and time set in file cap_restart,
#    and N is the global cubed-sphere grid resolution set in file
#    setCommonRunSettings.sh. A message is printed showing symlink and target
#    names. The program exits with an error message if the target file
#    is not found.
#
# USAGE: ./setRestartLink.sh
# 
####################################

rst_link_name=gchp_restart.nc4

# Get simulation start from cap_restart
if [ -f cap_restart ]; then
   start_str=$(sed 's/ /_/g' cap_restart)
else
   echo "ERROR: Unable to set ${rst_link_name} link because cap_restart does not exist! Create cap_restart containing simulation start date with format YYYYMMDD HHmmSS."
   exit	1
fi

# Set restart name, check that file exists, and set symlink
N=$(grep "CS_RES=" setCommonRunSettings.sh | cut -c 8- | xargs)
rst_target=./Restarts/GEOSChem.Restart.${start_str:0:13}z.c${N}.nc4
if [[ -f "${rst_target}" ]]; then
   ln -nsf ${rst_target} ${rst_link_name}
   echo "Restart symlink ${rst_link_name} set to ${rst_target}"
else
  echo "ERROR: Unable to set symlink ${rst_link_name} because file ${rst_target} does not exist! Create file or link with that name, or change start date in cap_restart and/or grid resolution in setCommonRunSettings.sh to match restart file that exists."
  exit 1
fi
