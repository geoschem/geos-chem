#!/bin/bash

#SBATCH -c 24
#SBATCH -N 1
#SBATCH -t 0-03:00
#SBATCH -p huce_intel,seas_compute,shared
#SBATCH --mem=30000
#SBATCH --mail-type=END

#############################################################################
### parallelTest.sh
###
### Performs a 1-hour GCClassic parallel test on Cannon within a 
### freshly-created run directory.  Uses 24 and 13 cores, which seems
### to reveal most parallelization issues.  Useful for when you only
### want to run a single parallel test instead of a whole suite of tests.
#############################################################################

#============================================================================
# %%% User-configurable settings %%%
#============================================================================
allCores=24                         # Run w/ same # cores as in #SBATCH -c
fewerCores=13                       # Run with fewer cores   
hhmm="0100"                         # Duration of run in hhmm
resultsLog="parallel.results.log"   # Logfile for results

#============================================================================
# %%% Bash functions %%%
#============================================================================

function get_subdir_name() {
    #========================================================================
    # Returns the subdirectory of OutputDir (or Restarts) where files
    # corresponding to a certain number of cores will be moved.
    #
    # 1st argument: Number of cores
    #========================================================================
    printf "${1}_cores"
    return 0
}


function find_restart_file() {

    #========================================================================
    # Returns the path of a GEOS-Chem restart file in a subdirectory
    #
    # 1st argument: Subdirectory (obtained from get_subdir_name)
    #========================================================================
    r=$(ls -1 -tr "Restarts/${1}" | grep "GEOS" | grep "${hhmm}" | head -1)
    printf "Restarts/${1}/${r}"
    return 0
}


function get_logfile_name() {
    #========================================================================
    # Returns the log file name corresponding to the number of cores
    #
    # 1st argument: Number of cores
    #========================================================================
    subdir=$(get_subdir_name "${1}")
    printf "parallel.${subdir}.log"
    return 0
}


function score_parallel_test() {
    #========================================================================
    # Determines if the parallelization test was successful by checking
    # that the restart files from both runs are bitwise identical.
    #
    # 1st argument: Folder where the 24 core restart file is located.
    # 2nd argument: Folder where the 13 core restart file is located.
    # 3rd argument: Name of the log file to print results
    #========================================================================

    # Arguments
    subDir1=$(get_subdir_name "${1}")
    subDir2=$(get_subdir_name "${2}")
    results="${3}"

    # Restart file names from both parallel test runs
    rstFile1=$(find_restart_file "${subDir1}")
    rstFile2=$(find_restart_file "${subDir2}")

    # Exit if eiher restart file does not exist
    [[ ! -f "${rstFile1}" ]] && return 1
    [[ ! -f "${rstFile2}" ]] && return 1

    # Remove the results log if it exists
    [[ -f "${results}" ]] && rm -f "${results}"

    # If the files are bitwise identical then the pqarallel test is successful
    diff "${rstFile1}" "${rstFile2}"
    if [[ $? -eq 0 ]]; then
	printf "\n\nParallel test result: PASS\n" >> "${results}"
	date                                      >> "${results}"
	return 0
    fi
    printf "\n\nParallel test result: FAIL\n"     >> "${results}"
    date                                          >> "${results}"
    return 1
}


function edit_config_files() {

    #========================================================================
    # Edits start times, end times, frequency, and duration in the various
    # run directory configuration files.
    #========================================================================

    # Replace ending time in geoschem_config.yml
    sed -i -e "s/20190801, ....../20190701, ${hhmm}00/" geoschem_config.yml
    sed -i -e "s/20190201, ....../20190101, ${hhmm}00/" geoschem_config.yml
    sed -i -e "s/20190201, ....../20190101, ${hhmm}00/" geoschem_config.yml
    sed -i -e "s/20130201, ....../20130101, ${hhmm}00/" geoschem_config.yml
    sed -i -e "s/20110201, ....../20110101, ${hhmm}00/" geoschem_config.yml

    # Replace freq & duration in HISTORY.rc
    sed -i -e "s/00000100 ....../00000000 ${hhmm}00/" HISTORY.rc

    # Change time cycle flag in HEMCO_Config.rc from EFYO to CYS,
    # to allow missing species to be set a default value.
    # Also make sure that 
    sed -i -e "s/EFYO/CYS/"            HEMCO_Config.rc  # GC_RESTART
    sed -i -e "s/EFY xyz 1/CYS xyz 1/" HEMCO_Config.rc  # GC_BCs
    
    # Make sure we get a HEMCO_diagnostics file in the output
    sed -i -e "s/DiagnFreq:                   Monthly/DiagnFreq:                   End/" HEMCO_Config.rc
    
    # Return success
    return 0
}


function move_files() {

    #========================================================================
    # Moves diagnostic or restart files to a subdirectory
    #
    # 1st argument: Directory name ("OutputDir" or "Restarts")
    #========================================================================

    # Arguments
    dirName="${1}"
    subDir="${2}"
    runLog="${3}"
    
    # Make the folder if it does not exist
    [[ ! -d "${dirName}/${subDir}" ]] && mkdir -p "${dirName}/${subDir}"

    # Restart files have end timestamp; diagnostics have start timestamp
    [[ "x${dirName}" == "xRestarts" ]] && stamp="${hhmm}" || stamp="0000"
	
    # Move ouptut
    for file in ${dirName}/*${stamp}*.nc*; do
	file="${file/${dirName}\/}"
	srcFile="${dirName}/${file}"
	trgFile="${dirName}/${subDir}/${file}"
	echo "Moving ${srcFile} to ${trgFile}" >> "${runLog}"
	mv "${srcFile}" "${trgFile}"
    done
}


function run_gcclassic() {

    #========================================================================
    # Runs a GEOS-Chem Classic simulation with the specified number of
    # OpenMP cores.  Copies diagnostics and restart files to the proper
    # subdirectory in order to facilitate comparison.
    #
    # 1st argument: Number of cores for OpenMP
    #========================================================================

    # Get the log file and subdirectory names from the # of cores
    runLog=$(get_logfile_name "${1}")
    subDir=$(get_subdir_name "${1}")
    
    # Remove any leftover files in the run dir
    ./cleanRunDir.sh --no-interactive  >> "${runLog}"

    # Run GEOS-Chem Classic
    export OMP_NUM_THREADS="${1}"
    export OMP_STACKSIZE="500m"
    echo "Now using ${OMP_NUM_THREADS} cores" >> "${runLog}"
    srun -c ${OMP_NUM_THREADS} ./gcclassic    >> "${runLog}"

    # Move files to the proper subdirectory for later comparison
    move_files "OutputDir" "${subDir}" "${runLog}"
    move_files "Restarts"  "${subDir}" "${runLog}"

    # Return success
    return 0
}


function main() {
    
    #========================================================================
    # Performs a parallelization test with the settings specified above.
    #========================================================================

    # Error check
    if [[ ! -f ./gcclassic ]]; then
	echo "Error...Could not find gcclassic executable!" >> $results
	exit 1
    fi

    # Perform the Parallel test
    edit_config_files                                  
    run_gcclassic       "${allCores}"
    run_gcclassic       "${fewerCores}" 
    score_parallel_test "${allCores}" "${fewerCores}" "${resultsLog}"
    return $?
}


#============================================================================
# Call the main program and return its status
#============================================================================
main $@
exit $?
