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
LINE="\n====================================================\n"
SED_INPUT_GEOS_1="s/20190801 000000/20190701 002000/"
SED_INPUT_GEOS_2="s/20190201 000000/20190101 002000/"
SED_HISTORY_RC="s/00000100 000000/00000000 002000/"


function absolute_path() {
    #========================================================================
    # Function to return the absolute path from a relative path
    #
    # 1st argument = relative path
    #========================================================================
    printf "$(cd "$(dirname "${1}")"; pwd -P)/$(basename "${1}")"
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
    printf "Removing leftover run directories and scripts:\n"
    for file in ${1}/*; do
	path=$(absolute_path ${file})
	printf " ... ${path}\n";
	rm -rf ${path}
    done
    unset path
}


function config_and_build() {
    #========================================================================
    # Configures and compiles GEOS-Chem for build_type=RELEASE
    #
    # 1st argument = Root folder for tests (w/ many rundirs etc)
    # 2nd argument = GEOS-Chem run directory name
    # 3rd argument = Log file where stderr and stdout will be redirected
    #========================================================================
    root=${1}
    rundir=${2}
    log=${3}

    # Echo information
    printf "${LINE}Now compiling ${rundir}${LINE}" >> ${log} 2>&1

    # Code configuration
    cd ${rundir}/build
    cmake ../CodeDir >> ${log} 2>&1
    if [[ $? -ne 0 ]]; then
	msg="Configuration of code in ${rundir} halted with error!\n"
	printf ${msg} >> ${log} 2>&1
	exit 100
    fi

    # Code compilation
    make -j install >> ${log} 2>&1
    if [[ $? -ne 0 ]]; then
	msg="Compilation of code in ${rundir} halted with error!\n"
	printf ${msg} >> ${log} 2>&1
	exit 110
    fi

    # Switch back to the root directory
    cd ${root}

    # Free variables
    unset log
    unset root
    unset rundir

    return 0
}


function run_gcclassic() {
    #========================================================================
    # Runs GEOS-Chem Classic
    #
    # 1st argument = root folder for tests (w/ many rundirs etc)
    # 2nd argument = GEOS-Chem run directory name
    # Configures and compiles GEOS-Chem for build_type=RELEASE
    #========================================================================
    root=${1}
    rundir=${2}
    log=${3}

    # Switch to the run directory
    cd ${root}/${rundir}

    # Echo information
    printf "${LINE}Now running GCClassic in ${rundir}${LINE}" >> ${log} 2>&1

    # Test if the executable exists
    if [[ !(-f ./gcclassic) ]]; then
	msg="The ./gcclassic executable was not found in ${rundir}!\n"
	printf ${msg} >> ${log} 2>&1
	exit 200
    fi

    # Run the code
    ./gcclassic >> ${log} 2>&1
    if [[ $? -ne 0 ]]; then
	msg="Program execution in ${rundir} halted with error!\n"
	printf ${msg} >> ${log} 2>&1
	exit 210
    fi

    # Switch back to the root directory for next iteration
    cd ${root}

    # Free variables
    unset log
    unset root
    unset rundir

    return 0
}
