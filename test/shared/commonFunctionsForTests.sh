#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: commonFunctionsForTests.sh
#
# !DESCRIPTION: Contains common functions that are used by the test scripts.
#\\
#\\
# !REVISION HISTORY:
#  03 Nov 2020 - R. Yantosca - Initial version
#  See the subsequent Git history with the gitk browser!
#EOP
#------------------------------------------------------------------------------
#BOC

# Global variables
FILL=$(printf '.%.0s' {1..44})
SEP_MAJOR=$(printf '=%.0s' {1..78})
SEP_MINOR=$(printf '\055%.0s' {1..78})
SED_INPUT_GEOS_1="s/End   YYYYMMDD, hhmmss  : 20190801 000000/End   YYYYMMDD, hhmmss  : 20190701 002000/"
SED_INPUT_GEOS_2="s/End   YYYYMMDD, hhmmss  : 20190201 000000/End   YYYYMMDD, hhmmss  : 20190101 002000/"
SED_INPUT_GEOS_3="s/Start YYYYMMDD, hhmmss  : 20160101 000000/Start YYYYMMDD, hhmmss  : 20190101 000000/"
SED_INPUT_GEOS_4="s/End   YYYYMMDD, hhmmss  : 20160201 000000/End   YYYYMMDD, hhmmss  : 20190101 010000/"
SED_INPUT_GEOS_5="s/End   YYYYMMDD, hhmmss  : 20160101 010000/End   YYYYMMDD, hhmmss  : 20190101 010000/"
SED_HISTORY_RC_1="s/00000100 000000/00000000 002000/"
SED_HISTORY_RC_2="s/7440000/010000/"
SED_RUN_CONFIG_1="s/20160101 000000/20190101 000000/"
SED_RUN_CONFIG_2="s/20160201 000000/20190101 010000/"
SED_RUN_CONFIG_3="s/20190201 000000/20190101 010000/"
SED_RUN_CONFIG_4="s/00000100 000000/00000000 010000/"
SED_RUN_CONFIG_5="s/7440000/010000/"
SED_RUN_CONFIG_6="s/1680000/010000/"
CMP_PASS_STR="Configure & Build......PASS"
CMP_FAIL_STR="Configure & Build......FAIL"
EXE_PASS_STR="Execute Simulation.....PASS"
EXE_FAIL_STR="Execute Simulation.....FAIL"
EXE_TBD_STR="Execute Simulation.....TBD"


function absolute_path() {
    #========================================================================
    # Returns the absolute path from a relative path
    #
    # 1st argument = relative path
    #========================================================================
    if [[ -d ${1} ]]; then
	absPath=$(readlink -f "${1}")   # If directory exists, use readlink
    else
	absPath="${1/#\~/$HOME}"        # Otherwise, replace ~ with $HOME
    fi
    echo "${absPath}"
}


function is_valid_rundir() {
    #========================================================================
    # Function that tests if a file is a valid GEOS-Chem run directory
    #
    # 1st argument: File or directory to be tested
    #========================================================================
    if [[ -d ${1} ]]; then
	if [[ -f ${1}/input.geos && -f ${1}/HEMCO_Config.rc ]]; then
	    echo "TRUE"
	    return
	fi
    fi
    echo "FALSE"
}


function is_gchp_rundir() {
    #========================================================================
    # Function that tests if a run directory is a GCHP run directory
    #
    # 1st argument: Directory to be tested
    #========================================================================
    expr=$(is_valid_rundir ${1})
    if [[ "x${expr}" == "xTRUE" ]]; then
	if [[ -f ${1}/CAP.rc ]]; then
	    echo "TRUE"
	    return
	fi
    fi
    echo "FALSE"
}


function count_rundirs() {
    #========================================================================
    # Returns the number of run directories in the root folder.
    #
    # 1st argument: Root directory containing several GC run directories
    #========================================================================
    root=${1}
    thisDir=$(pwd -P)

    # Count run directories
    cd ${root}
    declare -i x=0
    for rundir in *; do
	expr=$(is_valid_rundir ${rundir})
	if [[ "x${expr}" == "xTRUE" ]]; then
	    let x++
	fi
    done
    cd ${thisDir}

    # Cleanup and quit
    echo $x
}


function update_config_files() {
    #========================================================================
    # Function to replace text in rundir config files
    #
    # 1st argument = Integration tests root path
    # 2nd argument = GEOS-Chem run directory name
    #========================================================================
    root=${1}
    runDir=${2}

    # Replace text in input.geos
    sed -i -e "${SED_INPUT_GEOS_1}" ${root}/${runDir}/input.geos
    sed -i -e "${SED_INPUT_GEOS_2}" ${root}/${runDir}/input.geos
    sed -i -e "${SED_INPUT_GEOS_3}" ${root}/${runDir}/input.geos
    sed -i -e "${SED_INPUT_GEOS_4}" ${root}/${runDir}/input.geos
    sed -i -e "${SED_INPUT_GEOS_5}" ${root}/${runDir}/input.geos
    sed -i -e "${SED_HISTORY_RC_1}" ${root}/${runDir}/HISTORY.rc
    sed -i -e "${SED_HISTORY_RC_2}" ${root}/${runDir}/HISTORY.rc

    # For GCHP only
    expr=$(is_gchp_rundir "${root}/${runDir}")
    if [[ "x${expr}" == "xTRUE" ]]; then

	# Replace text in run.config.sh
	sed -i -e "${SED_RUN_CONFIG_1}" ${root}/${runDir}/runConfig.sh
	sed -i -e "${SED_RUN_CONFIG_2}" ${root}/${runDir}/runConfig.sh
	sed -i -e "${SED_RUN_CONFIG_3}" ${root}/${runDir}/runConfig.sh
	sed -i -e "${SED_RUN_CONFIG_4}" ${root}/${runDir}/runConfig.sh
	sed -i -e "${SED_RUN_CONFIG_5}" ${root}/${runDir}/runConfig.sh
	sed -i -e "${SED_RUN_CONFIG_6}" ${root}/${runDir}/runConfig.sh

	# Copy the run scripts
	cp ${root}/${runDir}/runScriptSamples/*.sh ${root}/${runDir}
    fi
}


function create_rundir() {
    #========================================================================
    # Calls createRunDir.sh with appropriate arguments
    #
    # 1st argument: Commands to forcefeed into createRunDir.sh
    # 2nd argument: Integration test root folder
    # 3rd argument: GEOS-Chem rundir to be created
    # 4th argument: Log file where stdout & stderr will be redirected
    #========================================================================
    cmd=${1}
    root=$(absolute_path ${2})
    runDir=${3}
    log=${4}

    # Create run dir for a short simulation
    printf " ... ${root}/${runDir}\n"
    printf ${cmd} | ./createRunDir.sh >> ${log} 2>&1
    update_config_files ${root} ${runDir}
}


function cleanup_files() {
    #========================================================================
    # Removes all files and directories in a root folder
    #
    # 1st argument = root folder for tests (w/ many rundirs etc)
    #========================================================================
    if [[ "x${1}" != "x" ]]; then

	# Give user a chance to avoid removing files
	printf "\nRemoving files and directories in ${1}:\n"
	printf "If this is OK, type 'yes to proceed or 'no' to quit:\n"
	read answer
	if [[ ! ${answer} =~ [Yy][Ee][Ss] ]]; then
	    printf "Exiting..."
	    return 1
	fi

	# Remove files and unlink links
	printf "Removing ...\n"
	for file in ${1}/*; do
	    if [[ -h ${file} ]]; then
		unlink ${file}
	    else
		path=$(absolute_path ${file})
		if [[ -e ${path} ]]; then
		    printf " ... ${path}\n";
		    rm -rf ${path}
		fi
		unset path
	    fi
	done
    fi
}


function print_to_log() {
    #========================================================================
    # Prints a message (single line) to a log file
    #
    # 1st argument: Message to be printed
    # 2nd argument: Log file where message will be sent
    #========================================================================
    printf "%s\n" "${1}" >> ${2} 2>&1
}


function gcclassic_exe_name() {
    #========================================================================
    # Returns the executable name given a directory name
    #
    # 1st argument: Directory name
    #========================================================================

    # Turn on case-insensitivity
    shopt -s nocasematch

    # Default executable name
    exeFileName="gcclassic"

    # Append a suffix to the executable file name for specific directories
    for suffix in apm bpch rrtmg tomas; do
	if [[ ${1} =~ ${suffix} ]]; then
	    exeFileName+=".${suffix}"
	    break
	fi
    done

    # Turn off case-insensitivity
    shopt -u nocasematch

    # Return
    echo "${exeFileName}"
}


function gcclassic_config_options() {
    #========================================================================
    # Returns the GCClassic configuration options given a directory name
    #
    # 1st argument: Directory name
    #========================================================================

    # Arguments
    dir=${1}
    baseOptions=${2}

    # Local variables
    exeFileName=$(gcclassic_exe_name ${dir})

    # Turn on case-insensitivity
    shopt -s nocasematch

    # Pick the proper build options
    if [[ ${dir} =~ "apm" ]]; then
	options="${baseOptions} -DAPM=y -DEXE_FILE_NAME=${exeFileName}"
    elif [[ ${dir} =~ "bpch" ]]; then
	options="${baseOptions} -DBPCH_DIAG=y -DEXE_FILE_NAME=${exeFileName}"
    elif [[ ${dir} =~ "rrtmg" ]]; then
	options="${baseOptions} -DRRTMG=y -DEXE_FILE_NAME=${exeFileName}"
    elif [[ ${dir} =~ "tomas" ]]; then
	options="${baseOptions} -DTOMAS=y -DTOMAS_BINS=15 -DBPCH_DIAG=y -DEXE_FILE_NAME=${exeFileName}"
    else
	options="${baseOptions}"
    fi

    # Turn off case-insensitivity
    shopt -u nocasematch

    # Return
    echo "${options}"
}


function gcclassic_compiletest_name() {
    #========================================================================
    # Returns the GCClassic configuration options given a directory name
    #
    # 1st argument: Directory name
    #========================================================================

    # Arguments
    dir=${1}

    # Turn on case-insensitivity
    shopt -s nocasematch

    # Pick the proper build options
    if [[ ${dir} =~ "apm" ]]; then
	result="GCClassic with APM"
    elif [[ ${dir} =~ "bpch" ]]; then
	result="GCClassic with BPCH diagnostics"
    elif [[ ${dir} =~ "rrtmg" ]]; then
	result="GCClassic with RRTMG"
    elif [[ ${dir} =~ "tomas" ]]; then
	result="GCClassic with TOMAS15"
    else
	result="GCClassic"
    fi

    # Turn off case-insensitivity
    shopt -u nocasematch

    # Return
    echo "${result}"
}


function build_gcclassic() {
    #========================================================================
    # Configures and compiles GEOS-Chem (Classic or GCHP) for
    # CMAKE_BUILD_TYPE=Debug.
    #
    # 1st argument = Root folder for tests (w/ many rundirs etc)
    # 2nd argument = GEOS_Chem Classic build directory
    # 3rd argument = Log file where stderr and stdout will be redirected
    # 4th argument = (OPTIONAL) Log file where results will be printed
    #========================================================================

    # Arguments
    root=$(absolute_path ${1})
    buildDir=${2}
    log=${3}
    results=${4}
    baseOptions=${5}

    # Local variables
    codeDir=${root}/CodeDir
    configOptions=$(gcclassic_config_options ${buildDir} "${baseOptions}")
    message=$(gcclassic_compiletest_name "${buildDir}")
    passMsg="$message${FILL:${#message}}.....${CMP_PASS_STR}"
    failMsg="$message${FILL:${#message}}.....${CMP_FAIL_STR}"

    #---------------------------------------
    # Code configuration
    #---------------------------------------
    cd ${buildDir}
    cmake ${codeDir} ${configOptions} >> ${log} 2>&1
    if [[ $? -ne 0 ]]; then
	if [[ "x${results}" != "x" ]]; then
	    print_to_log "${failMsg}" ${results}
	fi
	return 1
    fi

    #----------------------------------------
    # Code compilation and installation
    #----------------------------------------
    make -j install >> ${log} 2>&1
    if [[ $? -ne 0 ]]; then
	if [[ "x${results}" != "x" ]]; then
	    print_to_log "${failMsg}" ${results}
	fi
	return 1
    fi

    #----------------------------------------
    # Cleanup & quit
    #----------------------------------------

    # If we have gotten this far, the run passed,
    # so update the results log file accordingly
    print_to_log "${passMsg}" ${results}

    # Switch back to the root directory
    cd ${root}
}


function build_gchp() {
    #========================================================================
    # Configures and compiles GCHP.
    #
    # 1st argument = Root folder for tests (w/ many rundirs etc)
    # 2nd argument = GEOS_Chem Classic build directory
    # 3rd argument = Log file where stderr and stdout will be redirected
    # 4th argument = Log file where results will be printed
    # 5th argument = Compilation options
    #========================================================================

    # Arguments
    root=$(absolute_path ${1})
    buildDir=${2}
    log=${3}
    results=${4}
    options=${5}

    # Local variables
    codeDir=${root}/CodeDir
    message="GCHP"
    passMsg="$message${FILL:${#message}}.....${CMP_PASS_STR}"
    failMsg="$message${FILL:${#message}}.....${CMP_FAIL_STR}"

    #---------------------------------------
    # Code configuration
    #---------------------------------------
    cd ${buildDir}
    cmake ${codeDir} ${options} >> ${log} 2>&1
    if [[ $? -ne 0 ]]; then
	if [[ "x${results}" != "x" ]]; then
	    print_to_log "${failMsg}" ${results}
	fi
	return 1
    fi

    #----------------------------------------
    # Code compilation and installation
    #----------------------------------------
    make -j install >> ${log} 2>&1
    if [[ $? -ne 0 ]]; then
	if [[ "x${results}" != "x" ]]; then
	    print_to_log "${failMsg}" ${results}
	fi
	return 1
    fi

    #----------------------------------------
    # Cleanup & quit
    #----------------------------------------

    # If we have gotten this far, the run passed,
    # so update the results log file accordingly
    print_to_log "${passMsg}" ${results}

    # Switch back to the root directory
    cd ${root}
}


function submit_gchp_slurm_job() {
    #========================================================================
    # Submits a GCHP SLURM job script for the execution test phase.
    # Subsequent jobs are submitted as SLURM job dependencies.
    #
    # 1st argument = Root folder for tests (w/ many rundirs etc)
    # 2nd argument = GCHP run directory name
    # 3rd argument = Id from the previous SLURM job
    #========================================================================

    # Arguments
    root=${1}
    runDir=${2}
    jobId=${3}

    # Change to the run directory
    cd ${root}/${runDir}

    # Remove any leftover files in the run dir
    # (similar to createRunDir.sh, but do not ask user to confirm)
    rm -f cap_restart
    rm -f gcchem*
    rm -f *.rcx
    rm -f *~
    rm -f gchp.log
    rm -f HEMCO.log
    rm -f PET*.log
    rm -f multirun.log
    rm -f GC*.log
    rm -f log.dryrun*
    rm -f logfile.000000.out
    rm -f slurm-*
    rm -f 1
    rm -f EGRESS
    rm -f core.*
    rm -f ./OutputDir/*.nc*

    # Copy the executable here
    cp ${root}/build/bin/gchp .

    # Redirect the log file
    log="${root}/logs/execute.${runDir}.log"

    # Submit jobs (except the first) as a SLURM dependency
    if [[ "x${jobId}" == "xnone" ]]; then
	output=$(sbatch --export=ALL gchp.slurm.sh)
	output=($output)
	jobId=${output[3]}
    else
	output=$(sbatch --export=ALL --dependency=afterany:${jobId} gchp.slurm.sh)
	output=(${output})
	jobId=${output[3]}
    fi

    # Change to the root folder
    cd ${root}

    # Return the jobId for the next iteration
    echo "${jobId}"
}
