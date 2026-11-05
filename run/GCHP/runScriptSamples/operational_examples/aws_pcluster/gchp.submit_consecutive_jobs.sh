#!/bin/bash

# Use this bash script to submit multiple consecutive GCHP jobs (gchp.run) to SLURM.
# This allows breaking up long duration simulations into jobs of shorter duration.
# Set start date in cap_restart, set run duration in setCommonRunSettings.sh, and pass
# the number of runs as an argument. You must have gchp.run in your run directory.

# One argument: number of runs

numRuns=${1}
if [[ "x${numRuns}" == "x" ]]; then
    echo "ERROR: Specify number of runs as argument, e.g. ./gchp.submit_consecutive_jobs.sh 2"
    exit
fi

slog=gchp.submit.log

# Sanity check number of runs, start date, and duration per run

echo " " >> $slog
echo "Submitting ${numRuns} jobs at" `date` >> $slog
echo "Duration per job:  $(grep "Run_Duration=" setCommonRunSettings.sh | cut -c 14- | xargs)"  >> $slog
echo "Start of 1st job:  $(cat cap_restart)"  >> $slog

# Submit first job

msg=$(sbatch gchp.run)
echo "Executing command: sbatch gchp.run" >> $slog
echo $msg >> $slog
IFS=', ' IFS=', ' read -r -a msgarray <<< "$msg"
jobid=${msgarray[3]}

# Submit additional jobs
for i in $(seq 1 $((numRuns-1)));
do
    msg=$(sbatch --dependency=afterok:$jobid gchp.run)
    echo "Executing command: sbatch --dependency=afterok:$jobid gchp.run"  >> $slog
    echo $msg >> $slog
    IFS=', ' IFS=', ' read -r -a msgarray <<< "$msg"
    jobid=${msgarray[3]}
done
