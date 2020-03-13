#!/bin/bash

#SBATCH -n 30
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -t 0-0:30
#SBATCH -p huce_intel
#SBATCH --mem=MaxMemPerNode
#SBATCH --mail-type=ALL
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

# Resource request tips: 
#  (1) Use #SBATCH -n 6 to request 6 cores total across all nodes
#  (2) Use #SBATCH --exclusive to prevent other users from sharing nodes in use
#      (ability to use --exclusive may be disabled on some clusters)
#  (3) Use --mem=50G to request 50 Gigabytes per node, 
#      or use --mem=MaxMemPerNode to request all memory per node
#  (4) Performance is enhanced by requesting entire nodes and all memory
#      even if GCHP is run on fewer cores per node than available.
#
# See SLURM documentation for descriptions of all possible settings.
# Type 'man sbatch' at the command prompt to browse documentation.

#-------------------------------------------------------------------------
# Initial version: Lizzie Lundgren, 7/12/2018

# This run script is for use within gchp.multirun.sh which submits 
# multiple consecutive GCHP jobs to break up longer runs into smaller jobs of 
# shorter duration.
#
# NOTES:
#  1. This run script is used within gchp.multirun.sh and should not be 
#     used on its own. Use gchp.run instead if not doing multi-segmented runs.
#  2. All stdout is sent to slurm.jobid.out; each run will have a separate log.
#  3. Log file multirun.log contains info about all runs, including job ids
#     and cancellation information if a problem is encountered.
#  4. Unlike at the start of gchp.run for single segmented runs, cap_restart
#     is not deleted at the start of multi-segemented runs; the content is
#     sent to the multirun log and is used to determine if a GCHP run was
#     successful. If it does not exist at the end of the run, or if its
#     content did not change since the last run, then all jobs will be
#     cancelled and a notification about that will be sent to multirun.log.

# See SLURM documentation for descriptions of the above settings as well
# as many other settings you may use.

# Define GEOS-Chem log files, except stdout and stderr which are sent to
# slurm.jobid.out and slurm.jobid.err, as configured with the sbatch commands at
# the top of this script.
multirunlog="multirun.log"

# Make sure GCHP output restart file does not exist with the original name
# used by MAPL. Its present will cause GCHP run to fail. The output restart
# file is renamed after a successful run, but a file with the original name
# may persist if the previous run was unsuccessful.
if [[ -e gcchem_internal_checkpoint ]]; then
    rm gcchem_internal_checkpoint
fi

# Define function to cancel all jobs. This is done if cap_restart does
# not exist after a run, or if its date string was not updated
cancel_all_jobs()
{
   msg="Submitted batch job"
   msgs=($(grep -hr "${msg}" ./${multirunlog}))
   for jobid in "${msgs[@]}"
   do
      if [[ ! "${msg}" = *"${jobid}"* ]]; then
         scancel ${jobid}
      fi
   done
}

# Source your environment file. This requires first setting the gchp.env
# symbolic link using script setEnvironment in the run directory. 
# Be sure gchp.env points to the same file for both compilation and 
# running. You can copy or adapt sample environment files located in 
# ./envSamples subdirectory.
gchp_env=$(readlink -f gchp.env)
if [ ! -f ${gchp_env} ] 
then
   echo "ERROR: gchp.rc symbolic link is not set!"
   echo "Copy or adapt an environment file from the ./envSamples "
   echo "subdirectory prior to running. Then set the gchp.env "
   echo "symbolic link to point to it using ./setEnvironment."
   echo "Exiting."
   exit 1
fi
echo "WARNING: You are using environment settings in ${gchp_env}"
source ${gchp_env}

if [ -e cap_restart ]; then
   cap_rst=$(cat cap_restart)
   echo " "
   echo "Cap_restart prior to run: ${cap_rst}"
else
   cap_rst="none"
   echo " "
   echo "No cap_restart prior to run"
fi

source runConfig.sh > /dev/null

# Monthly diagnostic option automatically updates HISTORY.rc freq/duration
# with number of hours for the month of the run. Includes leap year handling.
# NOTE: This feature is enabled from runConfig.sh where Monthly_Diag is defined
if [ ${Monthly_Diag} = "1" ]; then
   if [ -e cap_restart ]; then    
      yr_str=${cap_rst:0:4}
      mo_str=$((10#${cap_rst:4:2}))
   else
      yr_str=${Start_Time:0:4}
      mo_str=$((10#${Start_Time:4:2}))
   fi
   feb_hrs=672
   if [ "$((yr_str%4))" = "0" ]; then
      if [ ! "$((yr_str%100))" = "0" ] || [ "$((yr_str%400))" = "0" ]; then
         feb_hrs=696   
      fi
      fi
   dpm=(744 $feb_hrs 744 720 744 720 744 744 720 744 720 744)
   hrs_str=${dpm[$((mo_str-1))]}
   sed -i "s|common_freq\=.*|common_freq\=\"${hrs_str}0000\"|" ./runConfig.sh
   sed -i "s|common_dur\=.*|common_dur\=\"${hrs_str}0000\"|"   ./runConfig.sh
fi

echo " "
echo "Settings from runConfig.sh: "
echo " "
source ./runConfig.sh

# Use the last run's output restart file as the input restart file if
# cap_restart exists (meaning there was a previous run). Note that unlike
# in single runs and the first run in the multirun series, the restart file
# is set here rather than from runConfig.sh.
if [ -e cap_restart ]; then
   restart_datetime=$(echo $(cat cap_restart) | sed 's/ /_/g')
   sed -i "s|GIGCchem_INTERNAL_RESTART_FILE\:.*|GIGCchem_INTERNAL_RESTART_FILE\:     +gcchem_internal_checkpoint.restart.${restart_datetime}.nc4|" ./GCHP.rc
fi

# Only run GCHP if runConfig.sh had no errors
if [[ $? == 0 ]]; then

   NX=$( grep NX GCHP.rc | awk '{print $2}' )
   NY=$( grep NY GCHP.rc | awk '{print $2}' )
   coreCount=$(( ${NX} * ${NY} ))
   # Note that $coreCount can be less than the requested cores in #SBATCH -n
   
   # Use SLURM to distribute tasks across nodes
   planeCount=$(( ${coreCount} / ${SLURM_NNODES} ))
   if [[ $(( ${coreCount} % ${SLURM_NNODES} )) > 0 ]]; then
   	${planeCount}=$(( ${planeCount} + 1 ))
   fi
   
   # Echo info from computational cores to log file for displaying results
   echo " "
   echo "Running GCHP:"
   echo "# of cores: ${coreCount}"
   echo "# of nodes: ${SLURM_NNODES}"
   echo "cores per nodes: ${planeCount}"
   
   # Echo start date
   echo '===> Run started at' `date`
   
   # Odyssey-specific setting to get around connection issues at high # cores
   export OMPI_MCL_btl=openib

   # Run the simulation
   # SLURM_NTASKS is #SBATCH -n and SLURM_NNODES is #SBATCH -N above
   time srun -n ${coreCount} -N ${SLURM_NNODES} -m plane=${planeCount} --mpi=pmix ./geos

   # Echo end date
   echo '===> Run ended at' `date`
fi

# Check that cap_restart exists. If it does not, cancel all remaining jobs
if [ -f cap_restart ]; then    
   new_cap_rst=$(cat cap_restart)
   echo " "
   echo "cap_restart after run: ${new_cap_rst}" | tee -a ${multirunlog}
   if [[ "${new_cap_rst}" = "${cap_rst}" || "${new_cap_rst}" = "" ]]; then
      echo " "
      echo "Error: cap_restart did not update to different date!" >> ${multirunlog}
      echo "Cancelling all jobs." >> ${multirunlog}
      cancel_all_jobs
   else
      # If all went well, rename the restart (checkpoint) file for clarity
      # and to enable reuse as a restart file. MAPL cannot read in a file
      # with the same name as the output checkpoint filename configured in
      # GCHP.rc.
      if [ -f cap_restart ]; then
         restart_datetime=$(echo $(cat cap_restart) | sed 's/ /_/g')
         mv gcchem_internal_checkpoint gcchem_internal_checkpoint.restart.${restart_datetime}.nc4
      fi
   fi
else
   echo " "
   echo "Error: cap_restart does not exist after GCHP run!" >> ${multirunlog}
   echo "Cancelling all jobs." >> ${multirunlog}
   cancel_all_jobs
fi

# Exit normally
exit 0

