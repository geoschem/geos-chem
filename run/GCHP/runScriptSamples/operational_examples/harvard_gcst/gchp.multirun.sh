#!/bin/bash

# Initial version: Lizzie Lundgren, 7/12/2018

# Use this bash script to submit multiple consecutive GCHP jobs to SLURM. 
# This allows breaking up long duration simulations into jobs of shorter 
# duration. This also enables outputting monthly diagnostics by updating
# the diagnostics frequency and duration to # of hours per month each run.
#
# Each job runs gchp.multirun.run once, with duration configured in 
# runConfig.sh. The first job submitted starts at the start date configured 
# in runConfig.sh. Subsequent jobs start where the last job left off, with 
# the new start date-time string stored in cap_restart. Since cap_restart is 
# over-written with every job, the start dates are sent to log file 
# cap_restart.log for later inspection. 
#
# Special Notes:
#
#  1. Configure the total number of runs via command line argument. This is 
#     equivalent to how many jobs will be submitted to SLURM. Make sure that 
#     the end date in runConfig.sh is sufficiently PAST the start date to 
#     accommodate all configured runs. It does not need to equal actual last
#     run end date.
#
#  2. This script uses a special run script for multi-run segments 
#     (gchp.multirun.run). Using the default run script instead (gchp.run) 
#     will not work for multi-segmented runs without updates. 
#
#  3. The run script submitted on loop in this shell script will send stdout to
#     to slurm-jobid.out. Log output is not over-written by subsequent runs. 
#
#  4. Because GCHP output diagnostic files contain the date, diagnostic output
#     files will also not be over-written by subsequent runs. 
#
#  5. Restart files will always be produced at the end of a run for use in 
#     the next run, but will not include the date in the filename. The run
#     script therefore renames them to include the date. If your run fails
#     mid-run then that renaming does not happen. In that case you may see 
#     a restart file called gcchem_internal_checkpoint. It will get deleted
#     upon rerun of the multi-run script. Properly renamed restart files
#     will not get deleted.
# 
#  6. You may use the archiveRun.sh script to archive ALL runs in the multi-run. However,
#     note that the config files saved will only reflect settings for the last run in the
#     series.

# Set argument to number of runs
Num_Runs=${1}
if [[ "x${Num_Runs}" == "x" ]]; then
    echo "ERROR: Specify number of runs as argument, e.g. ./gchp.multirun.sh 12"
    exit
fi

# Set multirun log filename (separate from GEOS-Chem log file gchp.log)
multirunlog="multirun.log"

rm -f gchp.log
rm -f $multirunlog

# Source runConfig.sh to update rc files and get variables used below.
source runConfig.sh > $multirunlog

# Only continue if runConfig.sh had no errors
if [[ $? == 0 ]]; then

   echo "Submitting    ${Num_Runs} jobs with duration '${Duration}'" | tee -a $multirunlog
   if [[ -e cap_restart ]]; then
       echo 'WARNING: cap_restart file exists. Starting simulation at:' | tee -a $multirunlog
       cat cap_restart | tee -a $multirunlog
   else
       echo "Start date:   ${Start_Time}" | tee -a $multirunlog
   fi
   echo "End date:     ${End_Time}" | tee -a $multirunlog
   echo "*** Check that end date is sufficiently past start date to span all intended runs ***"
   
   msg=$(sbatch gchp.multirun.run)
   echo $msg | tee -a $multirunlog
   IFS=', ' IFS=', ' read -r -a msgarray <<< "$msg"
   jobid=${msgarray[3]}
   
   for i in $(seq 1 $((Num_Runs-1))); 
   do
     msg=$(sbatch --dependency=afterok:$jobid gchp.multirun.run)
     echo $msg | tee -a $multirunlog
     IFS=', ' IFS=', ' read -r -a msgarray <<< "$msg"
     jobid=${msgarray[3]}
   done

else
   echo "Problem in sourcing runConfig.sh"
   cat $multirunlog
fi

exit 0

