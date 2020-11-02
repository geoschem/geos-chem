#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: createIntTests.sh
#
# !DESCRIPTION: Creates integration test run directories in a user-specified
#  root folder, and copies a run script there.
#\\
#\\
# !CALLING SEQUENCE:
#  ./createIntTests /path/to/integration/test/root/folder
#
# !REMARKS:
#  Right now we pass values to the existing ./createRunDir.sh,
#  but will implement a more elegant solution later.#
#
# !REVISION HISTORY:
#  09 Oct 2020 - R. Yantosca - Initial version
#  See the subsequent Git history with the gitk browser!
#EOP
#------------------------------------------------------------------------------
#BOC

#=============================================================================
# Check arguments
#=============================================================================
if [[ "x${1}" == "x" ]]; then
    echo "ERROR: The root-level directory for tests has not been specified!"
    exit
fi

#=============================================================================
# Local functions
#=============================================================================

function update_config_files() {
    # Function to replace text in rundir config files
    # 1st argument = Integration tests root path
    # 2nd argument = GEOS-Chem run directory name
    sed -i -e "${SED_INPUT_GEOS_1}" ${1}/${2}/input.geos
    sed -i -e "${SED_INPUT_GEOS_2}" ${1}/${2}/input.geos
    sed -i -e "${SED_HISTORY_RC}" ${1}/${2}/HISTORY.rc
}

function create_rundir() {
    # Function to call createRunDir with appropriate arguments
    # 1st argument = Commands to force-feed to createRunDir.sh
    # 2nd argument = Integration tests root directory path
    # 3rd argument = GEOS-Chem run directory name
    printf " ... ${2}/${3}\n"
    printf ${1} | ./createRunDir.sh >> ${LOG} 2>&1
    update_config_files ${2} ${3}
}

#=============================================================================
# Global variables
#=============================================================================
printf "====================================\n"
printf "Creating GEOS-Chem Integration Tests\n"
printf "====================================\n\n"

# Define variables
GCC_RUN_DIR="../../run/GCClassic"
THIS_DIR=`pwd -P`
ROOT=${1}

# sed string editor commands
SED_INPUT_GEOS_1="s/20190801 000000/20190701 002000/"
SED_INPUT_GEOS_2="s/20190201 000000/20190101 002000/"
SED_HISTORY_RC="s/00000100 000000/00000000 002000/"

# Log file
LOG=${ROOT}/logs/createIntTests.log

#=============================================================================
# Initial setup of integration test directory
#=============================================================================

# Remove leftover log file
rm -f ${LOG}

# Remove run directories in the test folder
printf "Removing leftover run directories:\n"
for DIR in ${ROOT}/*; do
    if [[ -d ${DIR} ]]; then
	echo " ... ${DIR}"
	rm -rf ${DIR}
    fi
done

# Copying the run script to the root test folder
printf "\nCopying run script to: ${ROOT}\n"
cp -f ./runIntTests.sh ${ROOT}

# Create log directory
if [[ !(-d ${ROOT}/logs) ]]; then
    printf "\nCreating log directory: ${ROOT}/logs\n"
    mkdir ${ROOT}/logs
fi

# Change to the directory where we will create the rundirs
cd $GCC_RUN_DIR
printf "\nCreating new run directories:\n"

#=============================================================================
# Create individual run directories: 4x5 - MERRA2 - 72L
#=============================================================================

DIR="merra2_4x5_benchmark"
create_rundir "1\n1\n2\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"    ${ROOT} ${DIR}

DIR="merra2_4x5_complexSOA"
create_rundir "1\n1\n3\n1\n1\n1\n1\n${ROOT}\n${DIR}\nn\n" ${ROOT} ${DIR}

DIR="merra2_4x5_complexSOA+SVPOA"
create_rundir "1\n1\n3\n2\n1\n1\n1\n${ROOT}\n${DIR}\nn\n" ${ROOT} ${DIR}

DIR="merra2_4x5_CH4"
create_rundir "3\n1\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"       ${ROOT} ${DIR}

# Not defined for 2x25 yet
#DIR="merra2_4x5_CO2"
#create_rundir "4\n1\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"       ${ROOT} ${DIR}

DIR="merra2_4x5_Hg"
create_rundir "5\n1\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"       ${ROOT} ${DIR}

DIR="merra2_4x5_POPs_BaP"
create_rundir "6\n1\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"       ${ROOT} ${DIR}

DIR="merra2_4x5_standard"
create_rundir "1\n1\n1\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"    ${ROOT} ${DIR}

DIR="merra2_4x5_standard+aciduptake"
create_rundir "1\n1\n5\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"    ${ROOT} ${DIR}

DIR="merra2_4x5_standard+APM"
create_rundir "1\n1\n7\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"    ${ROOT} ${DIR}

DIR="merra2_4x5_standard+marinePOA"
create_rundir "1\n1\n4\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"    ${ROOT} ${DIR}

DIR="merra2_4x5_standard+RRTMG"
create_rundir "1\n1\n8\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"    ${ROOT} ${DIR}

DIR="merra2_4x5_standard+TOMAS15"
create_rundir "1\n1\n6\n1\n1\n1\n1\n${ROOT}\n${DIR}\nn\n" ${ROOT} ${DIR}

DIR="merra2_4x5_tagCH4"
create_rundir "7\n1\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"       ${ROOT} ${DIR}

DIR="merra2_4x5_tagCO"
create_rundir "8\n1\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"       ${ROOT} ${DIR}

DIR="merra2_4x5_tagO3"
create_rundir "9\n1\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"       ${ROOT} ${DIR}

DIR="merra2_4x5_TransportTracers"
create_rundir "10\n1\n1\n1\n1\n${ROOT}\n${DIR}\nn\n"      ${ROOT} ${DIR}

#=============================================================================
# Cleanup and quit
#=============================================================================

# Switch back to the present directory
cd ${THIS_DIR}

# Free variables
unset GCC_RUN_DIR
unset GCHP_RUN_DIR
unset THIS_DIR
unset ROOT
unset SED_INPUT_GEOS_1
unset SED_INPUT_GEOS_2
unset SED_HISTORY_RC
unset LOG
unset DIR

#EOC
