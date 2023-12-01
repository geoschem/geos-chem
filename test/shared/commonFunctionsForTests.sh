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
#
# !REMARKS:
#  NOTE: Integration tests and parallelization tests run for 1 hour.
#  The exceptions are the nested-grid simulation tests, which run for 20
#  minutes due to the smaller dynamic timestep.
#EOP
#------------------------------------------------------------------------------
#BOC

# Global variables
FILL=$(printf '.%.0s' {1..47})
SEP_MAJOR=$(printf '=%.0s' {1..78})
SEP_MINOR=$(printf '\055%.0s' {1..78})
SED_CONFIG_1='s/end_date: \[20110201, 000000\]/end_date: \[20110101, 010000\]/'
SED_CONFIG_2='s/start_date: \[20160101, 000000\]/start_date: \[20190101, 000000\]/'
SED_CONFIG_3='s/end_date: \[20160201, 000000\]/end_date: \[20190101, 010000\]/'
SED_CONFIG_4='s/end_date: \[20160101, 000000\]/end_date: \[20190101, 010000\]/'
SED_CONFIG_5='s/end_date: \[20190201, 000000\]/end_date: \[20190101, 010000\]/'
SED_CONFIG_6='s/end_date: \[20190801, 000000\]/end_date: \[20190701, 010000\]/'
SED_CONFIG_N1='s/end_date: \[20190201, 000000\]/end_date: \[20190101, 002000\]/'
SED_CONFIG_N2='s/end_date: \[20190801, 000000\]/end_date: \[20190701, 002000\]/'
SED_HEMCO_CONF_1='s/GEOS_0.25x0.3125/GEOS_0.25x0.3125_NA/'
SED_HEMCO_CONF_2='s/GEOS_0.5x0.625/GEOS_0.5x0.625_NA/'
SED_HEMCO_CONF_3='s/DiagnFreq:                   Monthly/DiagnFreq:                   00000000 010000/'
SED_HEMCO_CONF_4='s/DiagnFreq:                   Monthly/DiagnFreq:                   00000000 002000/'
SED_HEMCO_CONF_N='s/\$RES.\$NC/\$RES.NA.\$NC/'
SED_HISTORY_RC_1='s/00000100 000000/00000000 010000/'
SED_HISTORY_RC_N='s/00000100 000000/00000000 002000/'
CMP_PASS_STR='Configure & Build.....PASS'
CMP_FAIL_STR='Configure & Build.....FAIL'
EXE_PASS_STR='Execute Simulation....PASS'
EXE_FAIL_STR='Execute Simulation....FAIL'
EXE_TBD_STR='Execute Simulation....TBD'
EXE_GCC_BUILD_LIST=("default" "apm"   "carbon"  "hg"       \
                    "luowd"   "rrtmg" "tomas15" "tomas40" )
EXE_GCHP_BUILD_LIST=("default" "carbon" "rrtmg")
END_hhmm_1h="0100z.nc4"
END_hhmm_20m="0020z.nc4"
PAR_TEST_SUFFIX="threads"
BIN_DIR="bin"
BUILD_DIR="build"
ENV_DIR="env"
LOGS_DIR="logs"
RUNDIRS_DIR="rundirs"
SCRIPTS_DIR="scripts"


function sed_ie() {
    #========================================================================
    # Replacement for `sed -i -e` that works on both MacOS and Linux
    #
    # 1st argument = regular expression
    # 2nd argument = file to be edited
    #========================================================================
    regex="${1}"
    file="${2}"
    if [[ "x$(uname -s)" == "xDarwin" ]]; then
        sed -i '' -e "${regex}" "${file}"          # MacOS/Darwin
    else
        sed -i -e "${regex}" "${file}"             # GNU/Linux
    fi
}


function sed_string() {
    #========================================================================
    # Returns a sed command string to replace a substring in a line of text.
    #
    # 1st argument = String of text
    # 2nd argument = Substring
    # 3rd argument = Replacement substring
    #========================================================================
    text=${1//\//\\\/}                # Replace '/' with '\/'
    text=${text//\*/\\\*}             # Replace '*' with '\*'
    text=${text//\./\\\.}             # Replace '.' with '\.'
    text=${text//\$/\\\$}             # Replace '$' with '\$'
    newText=${text/"${2}"/"${3}"}     # Repace substr w/ replacement str
    echo "s/${text}/${newText}/"  # Create sed command
}


function absolute_path() {
    #========================================================================
    # Returns the absolute path from a relative path
    #
    # 1st argument = relative path
    #========================================================================
    if [[ -d "${1}" ]]; then
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
    if [[ -d "${1}" ]]; then
        if [[ -f "${1}/geoschem_config.yml" && \
          -f "${1}/HEMCO_Config.rc" ]]; then
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
    expr=$(is_valid_rundir "${1}")
    if [[ "x${expr}" == "xTRUE" ]]; then
        if [[ -f "${1}/CAP.rc" ]]; then
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
    # 1st argument: Top-level folder where run directories are created
    #========================================================================
    rundirsDir="${1}"
    thisDir=$(pwd -P)

    # Count run directories
    cd "${rundirsDir}"
    declare -i x=0
    for runDir in *; do
        expr=$(is_valid_rundir "${runDir}")
        if [[ "x${expr}" == "xTRUE" ]]; then
            let x++
        fi
    done
    cd "${thisDir}"

    # Cleanup and quit
    echo $x
}


function update_gchp_config_files() {
    #========================================================================
    # Returns the number of run directories in the root folder.
    #
    # 1st argument: Root directory containing several GC run directories
    # 2nd argument: Integration test root directory
    #========================================================================
    runPath=$(absolute_path "${1}")           # GCHP run dir (abs path)

    # Edit config file setCommonRunSettings.sh
    file="${runPath}/setCommonRunSettings.sh"

    # 24 cores on 1 node
    sed_ie "s/TOTAL_CORES=.*/TOTAL_CORES=24/"                       "${file}"
    sed_ie "s/NUM_NODES=.*/NUM_NODES=1/"                            "${file}"
    sed_ie "s/NUM_CORES_PER_NODE=.*/NUM_CORES_PER_NODE=24/"         "${file}"

    # C24 grid resolution
    sed_ie "s/CS_RES=.*/CS_RES=24/"                                 "${file}"

    # 1-hr duration
    sed_ie "s/Run_Duration=\".*/Run_Duration=\"00000000 010000\"/"  "${file}"

    # 1-hr diagnostics
    sed_ie "s/AutoUpdate_Diagnostics=.*/AutoUpdate_Diagnostics=ON/" "${file}"
    sed_ie "s/Diag_Monthly=\"1\".*/Diag_Monthly=\"0\"/"             "${file}"
    sed_ie "s/Diag_Frequency=\".*/Diag_Frequency=\"010000\"/"       "${file}"
    sed_ie "s/Diag_Duration=\".*/Diag_Duration=\"010000\"/"         "${file}"

    # Do not require species in restart file
    sed_ie "s/Require_Species_in_Restart=.*/Require_Species_in_Restart=0/" "${file}"
}

function update_config_files() {
    #========================================================================
    # Function to replace text in rundir config files
    #
    # 1st argument = Integration tests root path
    # 2nd argument = GEOS-Chem run directory name
    #========================================================================
    runPath=$(absolute_path "${1}")     # GEOS-Chem rundir (abs path)

    #------------------------------------------------------------------------
    # Replace text in geoschem_config.yml
    #------------------------------------------------------------------------

    # For nested-grid fullchem runs, change simulation time to 20 minutes
    # in order to reduce the run time of the whole set of integration tests.
    if grep -q "05x0625" <<< "${runPath}"; then
        sed_ie "${SED_CONFIG_N1}" "${runPath}/geoschem_config.yml"
        sed_ie "${SED_CONFIG_N2}" "${runPath}/geoschem_config.yml"
    fi

    # Other text replacements
    sed_ie "${SED_CONFIG_1}" "${runPath}/geoschem_config.yml"
    sed_ie "${SED_CONFIG_2}" "${runPath}/geoschem_config.yml"
    sed_ie "${SED_CONFIG_3}" "${runPath}/geoschem_config.yml"
    sed_ie "${SED_CONFIG_4}" "${runPath}/geoschem_config.yml"
    sed_ie "${SED_CONFIG_5}" "${runPath}/geoschem_config.yml"
    sed_ie "${SED_CONFIG_6}" "${runPath}/geoschem_config.yml"

    #------------------------------------------------------------------------
    # Replace text in HEMCO_Config.rc
    #------------------------------------------------------------------------

    # For all nested-grid rundirs, add a NA into the entries for met fields
    # Also update the DiagnFreq for nested or global simulations
    if grep -q "05x0625" <<< "${runPath}"; then
        sed_ie "${SED_HEMCO_CONF_N}" "${runPath}/HEMCO_Config.rc"
        sed_ie "${SED_HEMCO_CONF_4}" "${runPath}/HEMCO_Config.rc"
    else
        sed_ie "${SED_HEMCO_CONF_3}" "${runPath}/HEMCO_Config.rc"
    fi

    # Other text replacements
    sed_ie "${SED_HEMCO_CONF_1}" "${runPath}/HEMCO_Config.rc"
    sed_ie "${SED_HEMCO_CONF_2}" "${runPath}/HEMCO_Config.rc"

    #------------------------------------------------------------------------
    # Replace text in HEMCO_Config.rc.gmao_metfields (GCClassic only)
    #------------------------------------------------------------------------

    if [[ -f "${runPath}/HEMCO_Config.rc.gmao_metfields" ]]; then
        # For all nested-grid rundirs, add a NA into the entries for met fields
        if grep -q "05x0625" <<< "${runPath}"; then
            sed_ie "${SED_HEMCO_CONF_N}" \
           "${runPath}/HEMCO_Config.rc.gmao_metfields"
        fi

        # Other text replacements
        sed_ie "${SED_HEMCO_CONF_1}" "${runPath}/HEMCO_Config.rc.gmao_metfields"
        sed_ie "${SED_HEMCO_CONF_2}" "${runPath}/HEMCO_Config.rc.gmao_metfields"
    fi

    #------------------------------------------------------------------------
    # Replace text in HISTORY.rc
    #------------------------------------------------------------------------

    # For nested-grid fullchem runs, change frequency and duration to 20 mins
    # in order to reduce the run time of the whole set of integration tests.
    if grep -q "05x0625" <<< "${runPath}"; then
        sed_ie "${SED_HISTORY_RC_N}" "${runPath}/HISTORY.rc"
    fi

    # Other text replacements
    sed_ie "${SED_HISTORY_RC_1}" "${runPath}/HISTORY.rc"
}


function create_rundir() {
    #========================================================================
    # Calls createRunDir.sh with appropriate arguments
    #
    # 1st argument: Commands to forcefeed into createRunDir.sh
    # 2nd argument: Log file where stdout & stderr will be redirected
    # 3rd argument: Integration test root directory
    #
    # NOTE: Run directories will be created with default names.
    #========================================================================
    cmd="${1}"                     # Command for createRunDir.sh
    log=$(absolute_path "${2}")    # Log file for output of createRunDir.sh
    logTmp="${log}.tmp"            # Temporary log file
    logTmp2="${log}.tmp.tmp"       # Temporary log file
    runPath=""                     # Path to run directory (abs path)

    # Create run dir for a short simulation
    printf "${cmd}" | ./createRunDir.sh > "${logTmp}" 2>&1

    # Get the absolute rundir path from the output of createRunDir.sh
    runPath=$(grep -E "Created" "${logTmp}")
    runPath=${runPath/Created /}
    printf " ... ${runPath}\n"

    # Now concatenate the temporary log into the log file
    # (Or create the log file if it doesn't yet exist)
    if [[ -f "${log}" ]]; then
        cat "${log}" "${logTmp}" >> "${logTmp2}"
        mv -f "${logTmp2}" "${log}"
        rm -f "${logTmp}"
    else
        mv "${logTmp}" "${log}"
    fi

    # Change start & end dates etc. in rundir configuration files
    update_config_files "${runPath}"

    # If this is a GCHP run directory, then also replace text in
    # GCHP-specific rundir configuration files etc.
    expr=$(is_gchp_rundir "${runPath}")
    if [[ "x${expr}" == "xTRUE" ]]; then
        update_gchp_config_files "${runPath}"
    fi

    # Remove temporary logs
    rm -rf "${logTmp} ${logTmp2}"
}


function cleanup_files() {
    #========================================================================
    # Removes all files and directories in a root folder
    #
    # 1st argument = root folder for tests (w/ many rundirs etc)
    #========================================================================
    itRoot="${1}"
    if [[ "x${itRoot}" != "x" ]]; then

    # Exit if directory is already empty
    if [[ ! $(ls -A "${itRoot}") ]]; then
        echo "${itRoot} is empty... nothing to clean up!"
        return 0
    fi

    # Give user a chance to avoid removing files
    printf "\nRemoving files and directories in ${itRoot}:\n"
    printf "If this is OK, type 'yes to proceed or 'no' to quit:\n"
    read answer
    if [[ ! "${answer}" =~ [Yy][Ee][Ss] ]]; then
        printf "Exiting...\n"
        exit 1
    fi

    # Remove files and unlink links
    printf "Removing ...\n"
    for entry in ${itRoot}/*; do
        printf " ... ${entry}\n";
        [[ -L ${entry} ]] && unlink ${entry}
        [[ -e ${entry} ]] && rm -rf ${entry}
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
    printf "%s\n" "${1}" >> "${2}" 2>&1
}


function exe_name() {
    #========================================================================
    # Returns the executable name given a directory name
    #
    # 1st argument: Model (GCClassic or GCHP)
    # 2nd argument: Build directory name
    #========================================================================
    model="${1}"
    buildDir="${2}"

    # Turn on case-insensitivity
    shopt -s nocasematch

    # Default executable name
    exeFileName="none"
    [[ "x${model}" == "xgcclassic" ]] && exeFileName="${model}"
    [[ "x${model}" == "xgchp"      ]] && exeFileName="${model}"

    # Append a suffix to the executable file name for specific directories
    if [[ "x${model}" == "xgcclassic" ]]; then
    for suffix in ${EXE_GCC_BUILD_LIST[@]}; do
        if [[ "${buildDir}" =~ "${suffix}" ]]; then
        exeFileName+=".${suffix}"
        break
        fi
    done
    elif [[ "x${model}" == "xgchp" ]]; then
    for suffix in ${EXE_GCHP_BUILD_LIST[@]}; do
        if [[ "${buildDir}" =~ "${suffix}" ]]; then
        exeFileName+=".${suffix}"
        break
        fi
    done
    fi

    # Turn off case-insensitivity
    shopt -u nocasematch

    # Return
    echo "${exeFileName}"
}


function config_options() {
    #========================================================================
    # Returns the configuration options given a directory name
    #
    # 1st argument: Model type ("gcclassic" or "gchp")
    # 2nd argument: Directory name
    # 3rd argument: Base compilation options
    #========================================================================

    # Arguments
    model="${1}"
    dir=$(basename "${2}")  # Only take last part of path
    baseOptions="${3}"

    # Local variables
    exeFileName=$(exe_name "${model}" "${dir}" )

    # Turn on case-insensitivity
    shopt -s nocasematch

    # Pick the proper build options
    if [[ ${dir} =~ "apm" ]]; then
        options="${baseOptions} -DAPM=y -DEXE_FILE_NAME=${exeFileName}"
    elif [[ ${dir} =~ "carbon" ]]; then
        options="${baseOptions} -DMECH=carbon -DEXE_FILE_NAME=${exeFileName}"
    elif [[ ${dir} =~ "hg" ]]; then
        options="${baseOptions} -DMECH=Hg -DEXE_FILE_NAME=${exeFileName}"
    elif [[ ${dir} =~ "luowd" ]]; then
        options="${baseOptions} -DLUO_WETDEP=y -DEXE_FILE_NAME=${exeFileName}"
    elif [[ ${dir} =~ "rrtmg" ]]; then
        options="${baseOptions} -DRRTMG=y -DEXE_FILE_NAME=${exeFileName}"
    elif [[ ${dir} =~ "tomas15" ]]; then
        options="${baseOptions} -DTOMAS=y -DTOMAS_BINS=15 -DBPCH_DIAG=y -DEXE_FILE_NAME=${exeFileName}"
    elif [[ ${dir} =~ "tomas40" ]]; then
        options="${baseOptions} -DTOMAS=y -DTOMAS_BINS=40 -DBPCH_DIAG=y -DEXE_FILE_NAME=${exeFileName}"
    else
        options="${baseOptions}"
    fi

    # Turn off case-insensitivity
    shopt -u nocasematch

    # Return
    echo "${options}"
}


function compiletest_name() {
    #========================================================================
    # Returns the GCClassic configuration options given a directory name
    #
    # 1st argument: Model name ('gcclassic' or 'gchp')
    # 2nd argument: Build directory name
    #========================================================================
    model="${1}"
    buildDir=$(basename "${2}")   # Only take last part of path

    # Display the proper model name
    [[ "x${model}" == "xgcclassic" ]] && displayName="GCClassic"
    [[ "x${model}" == "xgchp"      ]] && displayName="GCHP"

    # Turn on case-insensitivity
    shopt -s nocasematch

    # Pick the proper build options
    if [[ "${buildDir}" =~ "apm" ]]; then
        result="${displayName} with APM"
    elif [[ "${buildDir}" =~ "carbon" ]]; then
        result="${displayName} w/ carbon gases (as a KPP mechanism)"
    elif [[ "${buildDir}" =~ "luowd" ]]; then
        result="${displayName} with Luo et al wetdep"
    elif [[ "${buildDir}" =~ "hg" ]]; then
        result="${displayName} with Hg (as a KPP mechanism)"
    elif [[ "${buildDir}" =~ "rrtmg" ]]; then
        result="${displayName} with RRTMG"
    elif [[ "${buildDir}" =~ "tomas15" ]]; then
        result="${displayName} with TOMAS15"
    elif [[ "${buildDir}" =~ "tomas40" ]]; then
        result="${displayName} with TOMAS40"
    else
        result="${displayName}"
    fi

    # Turn off case-insensitivity
    shopt -u nocasematch

    # Return
    echo "${result}"
}


function build_model() {
    #========================================================================
    # Configures and compiles GEOS-Chem (Classic or GCHP) for
    # CMAKE_BUILD_TYPE=Debug.
    #
    # 1st argument = Model name ('gcclassic' or 'gchp')
    # 2nd argument = Root folder for tests (w/ many rundirs etc)
    # 3rd argument = CMake build directory
    # 4th argument = Base compilation options
    # 5th argument = Log file where stderr and stdout will be redirected
    # 6th argument = Log file where results will be printed
    #========================================================================

    # Arguments
    model="${1}"
    itRoot=$(absolute_path "${2}")
    thisBuildDir="${3}"
    baseOptions="${4}"
    log="${5}"
    results="${6}"

    # Stop with error if the model name is invalid
    if [[ "x${model}" != "xgcclassic" && "x${model}" != "xgchp" ]]; then
        echo "ERROR: '${model}' is an invalid model name!"
        echo "Exiting in 'build_model' (in 'commonFunctionsForTests.sh')"
        exit 1
    fi

    # Local variables
    codeDir="${itRoot}/CodeDir"
    configOptions=$(config_options "${model}" "${thisBuildDir}" "${baseOptions}")
    message=$(compiletest_name "${model}" "${thisBuildDir}")
    passMsg="$message${FILL:${#message}}.....${CMP_PASS_STR}"
    failMsg="$message${FILL:${#message}}.....${CMP_FAIL_STR}"

    # Change to the build directory
    cd "${thisBuildDir}"

    #---------------------------------------
    # Code configuration
    #---------------------------------------
    # Note: do not put quotes around configOptions, as bash will treat
    # it as a single string rather than individual options.
    cmake "${codeDir}" ${configOptions} >> "${log}" 2>&1
    if [[ $? -ne 0 ]]; then
        if [[ "x${results}" != "x" ]]; then
            print_to_log "${failMsg}" "${results}"
        fi
        return 1
    fi

    #----------------------------------------
    # Code compilation and installation
    #----------------------------------------
    make -j install >> "${log}" 2>&1
    if [[ $? -ne 0 ]]; then
        if [[ "x${results}" != "x" ]]; then
            print_to_log "${failMsg}" "${results}"
        fi
        return 1
    fi

    #----------------------------------------
    # Cleanup & quit
    #----------------------------------------

    # If we have gotten this far, the run passed,
    # so update the results log file accordingly
    print_to_log "${passMsg}" "${results}"
}


function rename_end_restart_file() {
    #========================================================================
    # Appends a suffix to the ending restart file in order to denote
    # the number of cores that were used.
    #
    # 1st argument: Number of OpenMP threads
    #========================================================================
    suffix="${1}${PAR_TEST_SUFFIX}"

    # Rename the ending restart file
    for r in Restarts/*.nc4; do
        if [[ "${r}" =~ "${END_hhmm_1h}" || "${r}" =~ "${END_hhmm_20m}" ]]; then
            newRstFile="${r}.${suffix}"
            mv -f "${r}" "${newRstFile}"
            [[ -f "${newRstFile}" ]] && return 0
        fi
    done
    return 1
}


function score_parallelization_test() {
    #========================================================================
    # Determines if the parallelization test was successful by checking
    # that the restart files from both runs are bitwise identical.
    #
    # 1st argument: Number of OpenMP threads used in 1st parallel test run
    # 2nd argument: Number of OpenMP threads used in 2nd parallel test run
    #========================================================================

    # Restart file names from both parallel test runs
    rstFile1=$(ls -1 Restarts/*.nc4* | grep "${1}${PAR_TEST_SUFFIX}")
    rstFile2=$(ls -1 Restarts/*.nc4* | grep "${2}${PAR_TEST_SUFFIX}")

    # Exit if eiher restart file does not exist
    [[ ! -f "$rstFile1" ]] && return 1
    [[ ! -f "$rstFile2" ]] && return 1

    # If the files are bitwise identical then the parallel test is successful
    diff "$rstFile1" "$rstFile2"
    return $?
}


function print_bootstrap_info_message() {
    #========================================================================
    # Prints a message to indicate if species bootstrapping is enabled.
    #
    # 1st argument: Enable ("yes") or disable ("no") bootstrapping
    #========================================================================
    if [[ "x${1}" == "xyes" ]]; then
        echo ""
        echo "%%%%% Species not found in the restart file will   %%%%%"
        echo "%%%%% be bootstrapped (i.e. set to default values) %%%%%"
    else
        echo ""
        echo "%%%%% Integration tests will fail unless all %%%%%"
        echo "%%%%% species are found in the restart file! %%%%%"
    fi
}


function change_time_cycle_flags() {
    #========================================================================
    # Changes the HEMCO time cycle flag for a given HEMCO container.
    #
    # 1st argument = HEMCO configuration file
    # 2nd argument = HEMCO container
    # 3rd argument = Old time cycle flag
    # 4th argument = New time cycle flag
    #========================================================================
    hcoCfg="${1}"
    hcoCont="${2}"
    oldFlag="${3}"
    newFlag="${4}"

    # Search for the container in HEMCO configuration file
    # The set -f prevents * from being expanded to a file listing
    set -f
    text=$(grep "${hcoCont}" "${hcoCfg}")
    unset -f

    # Replace the time cycle flag for the container
    if [[ "x${text}" != "x" ]]; then
	sedCmd=$(sed_string "${text}" "${oldFlag}" "${newFlag}")
	sed_ie "${sedCmd}" "${hcoCfg}"
    fi
}


function gcc_enable_or_disable_bootstrap() {
    #========================================================================
    # Edits HEMCO_Config.rc files to enable or disable "bootstrapping",
    # which is setting missing restart file variables to default values.
    #
    # 1st argument: Enable ("yes") or disable ("no") bootstrapping
    # 2nd argument: Directory containing integration test rundirs
    #========================================================================
    bootStrap="${1}"
    rundirsDir="${2}"

    print_bootstrap_info_message "${bootStrap}"

    # Loop over rundirs
    for runDir in ${rundirsDir}/*; do

        # Do the following if for only valid GEOS-Chem run dirs
        expr=$(is_valid_rundir "${runDir}")
        if [[ "x${expr}" == "xTRUE" ]]; then

            # Specify path HEMCO_Config.rc
            hcoCfg="${runDir}/HEMCO_Config.rc"

            if [[ "x${bootStrap}" == "xyes" ]]; then
		# Set missing species in restarts & BC files to defaults
		change_time_cycle_flags "${hcoCfg}" "SPC_ " "EFYO" "CYS"
		change_time_cycle_flags "${hcoCfg}" "SPC_ " "EY"   "CYS"
		change_time_cycle_flags "${hcoCfg}" "BC_ "  "EFY"  "CYS"
            else
		# Halt run if species are missing from restarts & BC files
		change_time_cycle_flags "${hcoCfg}" "SPC_ " "CYS"  "EFYO"
		change_time_cycle_flags "${hcoCfg}" "SPC_ " "EY"   "EFYO"
		change_time_cycle_flags "${hcoCfg}" "BC_ "  "CYS"  "EFY"
            fi
        fi
    done
}


function gchp_enable_or_disable_bootstrap() {
    #========================================================================
    # Edits HEMCO_Config.rc files to enable or disable "bootstrapping",
    # which is setting missing restart file variables to default values.
    #
    # 1st argument: Enable ("yes") or disable ("no") bootstrapping
    # 2nd argument: Directory containing integration test rundirs
    #========================================================================
    bootStrap="${1}"
    rundirsDir="${2}"

    print_bootstrap_info_message "${bootStrap}"

    # Loop over rundirs
    for runDir in ${rundirsDir}/*; do

        # Do the following if for only valid GEOS-Chem run dirs
        expr=$(is_gchp_rundir "${runDir}")
        if [[ "x${expr}" == "xTRUE" ]]; then

            # Specify path HEMCO_Config.rc
            script="${runDir}/setCommonRunSettings.sh"

            if [[ "x${bootStrap}" == "xyes" ]]; then
                # Set missing restart file variables to defaults
                sed_ie "s/Require_Species_in_Restart=./Require_Species_in_Restart=0/" "${script}"
            else
                # Don't set missing restart file variables to defaults
                sed_ie "s/Require_Species_in_Restart=./Require_Species_in_Restart=1/" "${script}"
            fi
        fi
    done
}
