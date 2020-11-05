
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
SED_INPUT_GEOS_1="s/20190801 000000/20190701 002000/"
SED_INPUT_GEOS_2="s/20190201 000000/20190101 002000/"
SED_HISTORY_RC="s/00000100 000000/00000000 002000/"
CMP_PASS_STR="Configure & Build......PASS"
CMP_FAIL_STR="Configure & Build......FAIL"
EXE_PASS_STR="Execute Simulation.....PASS"
EXE_FAIL_STR="Execute Simulation.....FAIL"


function absolute_path() {
    #========================================================================
    # Function to return the absolute path from a relative path
    #
    # 1st argument = relative path
    #========================================================================
    printf "$(cd "$(dirname "${1}")"; pwd -P)/$(basename "${1}")"
}


function count_rundirs() {
    #========================================================================
    # Returns the number of run directories in the root folder.
    #
    # 1st argument: Root directory containing several GC run directories
    #========================================================================
    root=${1}
    this_dir=$(pwd -P)

    # Count run directories
    cd ${root}
    declare -i x=0
    for rundir in *; do
	if [[ -d ${rundir} && "x${rundir}" != "xlogs" ]]; then
	    let x++
	fi
    done
    cd ${this_dir}

    # Cleanup and quit
    echo $x
    unset count
    unset root
    unset this_dir
}


function update_config_files() {
    #========================================================================
    # Function to replace text in rundir config files
    #
    # 1st argument = Integration tests root path
    # 2nd argument = GEOS-Chem run directory name
    #========================================================================
    root=${1}
    rundir=${2}

    # Replace text
    sed -i -e "${SED_INPUT_GEOS_1}" ${root}/${rundir}/input.geos
    sed -i -e "${SED_INPUT_GEOS_2}" ${root}/${rundir}/input.geos
    sed -i -e "${SED_HISTORY_RC}" ${root}/${rundir}/HISTORY.rc

    # Free variables
    unset root
    unset rundir
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
    rundir=${3}
    log=${4}

    # Create run dir for a short simulation
    printf " ... ${root}/${rundir}\n"
    printf ${cmd} | ./createRunDir.sh >> ${log} 2>&1
    update_config_files ${root} ${rundir}

    # Free variables
    unset cmd
    unset log
    unset root
    unset rundir
}


function cleanup_files() {
    #========================================================================
    # Removes all files and directories in a root folder
    #
    # 1st argument = root folder for tests (w/ many rundirs etc)
    #========================================================================
    if [[ "x${1}" != "x" ]]; then
	printf "Removing leftover run directories and scripts:\n"
	for file in ${1}/*; do
	    path=$(absolute_path ${file})
	    printf " ... ${path}\n";
	    rm -rf ${path}
	done
	unset path
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


function count_matches_in_log() {
    #========================================================================
    # Computes the number of times a string is found in a log file.
    #
    # 1st argument: Search string
    # 2nd argument: Log file
    #========================================================================
    matches=$(grep "${1}" ${2} | wc -l)
    echo "${matches}"
}


function config_and_build() {
    #========================================================================
    # Configures and compiles GEOS-Chem for build_type=RELEASE
    #
    # 1st argument = Root folder for tests (w/ many rundirs etc)
    # 2nd argument = GEOS-Chem run directory name
    # 3rd argument = Log file where stderr and stdout will be redirected
    # 4th argument = (OPTIONAL) Log file where results will be printed
    #========================================================================

    # Local variables
    root=${1}
    rundir=${2}
    log=${3}
    results=${4}
    passMsg="$rundir${FILL:${#rundir}}.....${CMP_PASS_STR}"
    failMsg="$rundir${FILL:${#rundir}}.....${CMP_FAIL_STR}"

    #---------------------
    # Code configuration
    #---------------------
    cd ${rundir}/build
    cmake ../CodeDir >> ${log} 2>&1
    if [[ $? -ne 0 ]]; then
	if [[ "x${results}" != "x" ]]; then
	    print_to_log "${failMsg}" ${results}
	fi
	return 1
    fi

    #---------------------
    # Code compilation
    #---------------------
    make -j install >> ${log} 2>&1
    if [[ $? -ne 0 ]]; then
	if [[ "x${results}" != "x" ]]; then
	    print_to_log "${failMsg}" ${results}
	fi
	return 1
    fi

    #---------------------
    # Cleanup & quit
    #---------------------

    # If we have gotten this far, the run passed,
    # so update the results log file accordingly
    print_to_log "${passMsg}" ${results}

    # Switch back to the root directory
    cd ${root}

    # Free variables
    unset failMsg
    unset log
    unset passMsg
    unset results
    unset root
    unset rundir
}


function run_gcclassic() {
    #========================================================================
    # Runs GEOS-Chem Classic
    #
    # 1st argument = root folder for tests (w/ many rundirs etc)
    # 2nd argument = GEOS-Chem run directory name
    # 3rd argument = Log file where stderr and stdout will be redirected
    # 4th argument = (OPTIONAL) Log file where results will be printed
    #========================================================================

    # Local variables
    root=${1}
    rundir=${2}
    log=${3}
    results=${4}
    passMsg="$rundir${FILL:${#rundir}}.....${EXE_PASS_STR}"
    failMsg="$rundir${FILL:${#rundir}}.....${EXE_FAIL_STR}"

    # Switch to the run directory
    cd ${root}/${rundir}

    #--------------------------------
    # Test if the executable exists
    #--------------------------------
    if [[ !(-f ./gcclassic) ]]; then
	if [[ "x${results}" != "x" ]]; then
	    print_to_logs "${failMsg}" ${results}
	fi
	return 1
    fi

    # Run the code
    ./gcclassic >> ${log} 2>&1
    if [[ $? -ne 0 ]]; then
	if [[ "x${results}" != "x" ]]; then
	    print_to_logs "${failMsg}" ${results}
	fi
	return 1
    fi

    #--------------------------------
    # Cleanup & quit
    #--------------------------------

    # If we have gotten this far, the run passed,
    # so update the results log file accordingly
    print_to_log "${passMsg}" ${results}

    # Switch back to the root directory for next iteration
    cd ${root}

    # Free variables
    unset log
    unset root
    unset rundir
    unset results
}
