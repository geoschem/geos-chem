#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: setupForRestarts.sh
#
# !DESCRIPTION: Contains bash functions for copying GEOS-Chem Classic
#  restart files from a remote location to the run directory.  The restart
#  file paths are read from the relevant download_data.yml file so that
#  these paths can be updated in a single location.
#\\
#\\
# !REVISION HISTORY:
#  25 Nov 2023 - Initial version - R. Yantosca
#  See the Git history for additional updates
#------------------------------------------------------------------------------
#BOC

function parseYaml() {
    #========================================================================
    # Portable bash YAML parser by Stefan Farestam
    # See: stackoverflow.com/questions/5014632/how-can-i-parse-a-yaml-file-from-a-linux-shell-script
    #
    # 1st argument: YAML file to be parsed
    #
    # Usage: list=$(parse_yaml myfile.yaml)
    #========================================================================
    local prefix=$2
    local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
    sed -ne "s|^\($s\):|\1|" \
         -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
         -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
    awk -F$fs '{
       indent = length($1)/2;
       vname[indent] = $2;
       for (i in vname) {if (i > indent) {delete vname[i]}}
       if (length($3) > 0) {
          vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
          printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
       }
    }'
    return $?
}


function removeQuotes() {
    #========================================================================
    # Removes double quotes from a string
    #
    # 1st argument: String to be parsed
    #========================================================================
    echo ${1//\"/}
    return $?
}


function join() {
    #========================================================================
    # Joins a directory with a file name and prints the result
    #
    # 1st argument: Directory
    # 2nd argument: File
    #========================================================================
    dir=$(removeQuotes "${1}")
    file=$(removeQuotes "${2}")
    path="${dir}/${file}"
    echo $path
    return $?
}


function setRestartEnvVar() {
    #========================================================================
    # Exports environment variables for restart files that can
    # be used in createRunDir.sh and setupForRestarts.sh
    #
    # 1st argument: One line of output from parseYaml
    #========================================================================
    if [[ "x${1}" =~ "RUNDIR_restarts_" ]]; then
        export "${1}"
    fi
    return $?
}

function unsetRestartEnvVar() {
    #========================================================================
    # Unsets environment variables defined by setRestartEnvVar
    #
    # 1st argument: One line of output from parseYaml
    #========================================================================
    if [[ "x${1}" =~ "RUNDIR_restarts_" ]]; then
        array=(${1//\=/ })
        unset "${array[0]}"
    fi
    return $?
}


function getRemoteRoot() {
    #========================================================================
    # Returns the remote root directory for restart files
    # (i.e. all restart files are stored in subdirs of this dir)
    #
    # 1st argument: Are we on the AWS cloud?
    #========================================================================

    # Remote root path on AWS cloud
    if [[ "x${1}" != "x" ]]; then
        echo $(removeQuotes "s3://gcgrid/${RUNDIR_restarts_root}")
	return $?
    fi

    # Remote root path on a local cluster
    echo $(removeQuotes "${GC_DATA_ROOT}/${RUNDIR_restarts_root}")
    return $?
}


function getS3CopyCmd() {
    #========================================================================
    # Returns the AWS "s3 cp" command to be used
    #
    # 1st argument: Are we on the AWS cloud?
    #========================================================================

    # We are on the AWS cloud
    if [[ "x${1}" != "x" ]]; then
        echo "aws s3 cp --request-payer requester "
        return $?
    fi

    # We are not on the AWS cloud
    echo ""
    return $?
}


function getAerosolLocal() {
    #========================================================================
    # Returns local restart file path for the aerosol simulation
    #
    # 1st argument: Local restart file directory
    #========================================================================
    echo $(join "${1}" "${RUNDIR_restarts_aerosol_local}")
    return $?
}


function getAerosolRemote() {
    #========================================================================
    # Returns remote restart file path for the aerosol simulation
    #
    # 1st argument: Remote restart file directory
    #========================================================================
    echo $(join "${1}" "${RUNDIR_restarts_aerosol_remote}")
    return $?
}


function getCarbonLocal() {
    #========================================================================
    # Returns local restart file path for the carbon simulation
    #
    # 1st argument: Local restart file directory
    #========================================================================
    echo $(join "${1}" "${RUNDIR_restarts_carbon_local}")
    return $?
}


function getCarbonRemote() {
    #========================================================================
    # Returns remote restart file path for the carbon simulation
    #
    # 1st argument: Remote restart file directory
    #========================================================================
    echo $(join "${1}" "${RUNDIR_restarts_carbon_remote}")
    return $?
}


function getCH4Local() {
    #========================================================================
    # Returns local restart file path for the CH4 simulation
    #
    # 1st argument: Local restart file directory
    #========================================================================
    echo $(join "${1}" "${RUNDIR_restarts_ch4_local}")
    return $?
}


function getCH4Remote() {
    #========================================================================
    # Returns remote restart file path for the CH4 simulation
    #
    # 1st argument: Remote restart file directory
    #========================================================================
    echo $(join "${1}" "${RUNDIR_restarts_ch4_remote}")
    return $?
}


function getCO2Local() {
    #========================================================================
    # Returns local restart file path for the CO2 simulation
    #
    # 1st argument: Local restart file directory
    #========================================================================
    echo $(join "${1}" "${RUNDIR_restarts_co2_local}")
    return $?
}


function getCO2Remote() {
    #========================================================================
    # Returns remote restart file path for the CO2 simulation
    #
    # 1st argument: Remote restart file directory
    #========================================================================
    echo $(join "${1}" "${RUNDIR_restarts_co2_remote}")
    return $?
}


function getFullchemLocal() {
    #========================================================================
    # Returns local restart file path for fullchem simulations
    #
    # 1st argument: sim_extra_option
    # 2nd argument: Local restart file directory
    #========================================================================

    # TOMAS15
    if [[ "x${1}" == "xTOMAS15" ]]; then
        echo $(join "${2}" "${RUNDIR_restarts_tomas15_local}")
        return $?
    fi

    # TOMAS40
    if [[ "x${1}" == "xTOMAS40" ]]; then
        echo $(join "${2}" "${RUNDIR_restarts_tomas40_local}")
        return $?
    fi

    # Default fullchem
    echo $(join "${2}" "${RUNDIR_restarts_fullchem_local}")
    return $?
}


function getFullchemRemote() {
    #========================================================================
    # Returns local restart file path for fullchem simulations
    #
    # 1st argument: sim_extra_option
    # 2nd argument: Local restart file directory
    #========================================================================

    # TOMAS15
    if [[ "x${1}" == "xTOMAS15" ]]; then
       echo $(join "${2}" "${RUNDIR_restarts_tomas15_remote}")
       return $?
    fi

    # TOMAS40
    if [[ "x${1}" == "xTOMAS40" ]]; then
        echo $(join "${2}" "${RUNDIR_restarts_tomas40_remote}")
        return $?
    fi

    # Default fullchem
    echo $(join "${2}" "${RUNDIR_restarts_fullchem_remote}")
    return $?
}


function getMercuryLocal() {
    #========================================================================
    # Returns local restart file path for the mercury simulation
    #
    # 1st argument: Local restart file directory
    #========================================================================
    echo $(join "${1}" "${RUNDIR_restarts_mercury_local}")
    return $?
}


function getMercuryRemote() {
    #========================================================================
    # Returns remote restart file path for the mercury simulation
    #
    # 1st argument: Remote restart file directory
    #========================================================================
    echo $(join "${1}" "${RUNDIR_restarts_mercury_remote}")
    return $?
}


function getMetalsLocal() {
    #========================================================================
    # Returns local restart file path for the metals simulation
    #
    # 1st argument: Local restart file directory
    #========================================================================
    echo $(join "${1}" "${RUNDIR_restarts_metals_local}")
    return $?
}

function getMetalsRemote() {
    #========================================================================
    # Returns remote restart file path for the metals simulation
    #
    # 1st argument: Remote restart file directory
    #========================================================================
    echo $(join "${1}" "${RUNDIR_restarts_metals_remote}")
    return $?
}


function getPOPsLocal() {
    #========================================================================
    # Returns local restart file path for the POPs simulation
    #
    # 1st argument: Local restart file directory
    #========================================================================
    echo $(join "${1}" "${RUNDIR_restarts_pops_local}")
    return $?
}


function getPOPsRemote() {
    #========================================================================
    # Returns remote restart file path for the TagCH4 simulation
    #
    # 1st argument: Remote restart file directory
    #========================================================================
    echo $(join "${1}" "${RUNDIR_restarts_pops_remote}")
    return $?
}


function getTagCH4Local() {
    #========================================================================
    # Returns local restart file path for the TagCH4 simulation
    #
    # 1st argument: Local restart file directory
    #========================================================================
    echo $(join "${1}" "${RUNDIR_restarts_tagch4_local}")
    return $?
}


function getTagCH4Remote() {
    #========================================================================
    # Returns remote restart file path for the TagCH4 simulation
    #
    # 1st argument: Remote restart file directory
    #========================================================================
    echo $(join "${1}" "${RUNDIR_restarts_tagch4_remote}")
    return $?
}


function getTagCOLocal() {
    #========================================================================
    # Returns local restart file path for the TagCO simulation
    #
    # 1st argument: Local restart file directory
    #========================================================================
    echo $(join "${1}" "${RUNDIR_restarts_tagco_local}")
    return $?
}


function getTagCORemote() {
    #========================================================================
    # Returns remote restart file path for the TagCO simulation
    #
    # 1st argument: Remote restart file directory
    #========================================================================
    echo $(join "${1}" "${RUNDIR_restarts_tagco_remote}")
    return $?
}


function getTagO3Local() {
    #========================================================================
    # Returns local restart file path for the TagO3 simulation
    #
    # 1st argument: Local restart file directory
    #========================================================================
    echo $(join "${1}" "${RUNDIR_restarts_tago3_local}")
    return $?
}


function getTagO3Remote() {
    #========================================================================
    # Returns remote restart file path for the TagO3 simulation
    #
    # 1st argument: Remote restart file directory
    #========================================================================
    echo $(join "${1}" "${RUNDIR_restarts_tago3_remote}")
    return $?
}


function getTracersLocal() {
    #========================================================================
    # Returns local restart file path for the TransportTracers simulation
    #
    # 1st argument: Local restart file directory
    #========================================================================
    echo $(join "${1}" "${RUNDIR_restarts_transporttracers_local}")
    return $?
}


function getTracersRemote() {
    #========================================================================
    # Returns remote restart file path for the TransportTracers simulation
    #
    # 1st argument: Remote restart file directory
    #========================================================================
    echo $(join "${1}" "${RUNDIR_restarts_transporttracers_remote}")
    return $?
}


function copyRemoteToLocal() {
    #========================================================================
    # Copy the remote restart file to the local run directory.
    #
    remote_rst="${1}"   # 1st argument: Remote restart file path
    local_rst="${2}"    # 2nd argument: Local restart file path
    is_aws="${3}"       # 3rd argument: Are we on the AWS cloud?
    s3_cp="${4}"        # 4th argument: AWS "s3 cp" command to execute
    #========================================================================

    # If we are on AWS, copy the remote restart file from s3://gcgrid
    if [[ "x${is_aws}" != "x" ]]; then
        ${s3_cp} "${remote_rst}" "${local_rst}"
        return $?
    fi

    # Otherwise copy the remote restart file from the remote root directory
    if [[ -f "${remote_rst}" ]]; then
        cp "${remote_rst}" "${local_rst}"
        return $?
    fi

    # Otherwise note that the sample restart file was not found
    printf "\n  -- The following sample restart provided for this simulation was not found:"
    printf "\n     ${remote_rst}"
    printf "\n     You will need to provide this initial restart file or disable"
    printf "\n     GC_RESTARTS in HEMCO_Config.rc to initialize your simulation"
    printf "\n     with default background species concentrations.\n"
    return $?
}


function setEFYOtoEYinHemcoConfig() {
    #========================================================================
    # Sample restarts for several simulations do not contain all species.
    # For those simulations, print a warning and change the time
    # cycle option in HEMCO config so that we do not force an error
    # if not found (i.e. EFYO --> EY)
    #
    sim_name="${1}"           # 1st argument: Simulation name
    sim_extra_option="${2}"   # 2nd argument: Simulation extra option
    #========================================================================

    if [[ "x${sim_extra_option}" == "xaciduptake"       ||
          "x${sim_extra_option}" == "xmarinePOA"        ||
          "x${sim_extra_option}" == "xcomplexSOA_SVPOA" ||
          "x${sim_extra_option}" == "xAPM"              ||
          "x${sim_name}"         == "xPOPs"             ||
          "x${sim_name}"         == "xtagCH4"           ||
          "x${sim_name}"         == "xtagO3"        ]]; then
        old="SpeciesRst_?ALL?    \$YYYY/\$MM/\$DD/\$HH EFYO"
        new="SpeciesRst_?ALL?    \$YYYY/\$MM/\$DD/\$HH EY  "
        sed_ie "s|${old}|${new}|" HEMCO_Config.rc

        printf "\n  -- The sample restart provided for this simulation may not"
        printf "\n     contain all species defined in this simulation. Missing"
        printf "\n     species will be assigned default background concentrations."
        printf "\n     Check your GEOS-Chem log file for details. As always, it"
        printf "\n     is recommended that you spin up your simulation to ensure"
        printf "\n     proper initial conditions.\n"
    fi
    return $?
}


function copyRestartToRunDir() {
    #========================================================================
    # Copies a restart file to a run directory for a given
    # GEOS-Chem Classic simulation.
    #
    sim_name="${1}"           # 1st argument: Simulation name
    sim_extra_option="${2}"   # 2nd argument: Extra option for simulation
    rst_root="${3}"           # 3rd argument: Remote restart file dir
    loc_root="${4}"           # 4th argument: Local restart file dir
    #========================================================================

    # Determine the remote and local restart file paths for each simulation
    if [[ "x${sim_name}" == "xaerosol" ]]; then
        remote_rst=$(getAerosolRemote "${rst_root}")
        local_rst=$(getAerosolLocal "${loc_root}")

    elif [[ "x${sim_name}" == "xcarbon" ]]; then
        remote_rst=$(getCarbonRemote "${rst_root}")
	local_rst=$(getCarbonLocal "${loc_root}")

    elif [[ "x${sim_name}" == "xCH4" ]]; then
        remote_rst=$(getCH4Remote "${rst_root}")
        local_rst=$(getCH4Local "${loc_root}")

    elif [[ "x${sim_name}" == "xCO2" ]]; then
        remote_rst=$(getCO2Remote "${rst_root}")
        local_rst=$(getCO2Local "${loc_root}")

    elif [[ "x${sim_name}" == "xfullchem" ]]; then
        # NOTE: Also handles TOMAS40 and TOMAS15
        remote_rst=$(getFullchemRemote "${sim_extra_option}" "${rst_root}")
        local_rst=$(getFullchemLocal "${sim_extra_option}" "${loc_root}")

    elif [[ "x${sim_name}" == "xmetals" ]]; then
        remote_rst=$(getMetalsRemote "${rst_root}")
        local_rst=$(getMetalsLocal "${loc_root}")

    elif [[ "x${sim_name}" == "xHg" ]]; then
        remote_rst=$(getMercuryRemote "${rst_root}")
        local_rst=$(getMercuryLocal "${loc_root}")

    elif [[ "x${sim_name}" == "xPOPs" ]]; then
        remote_rst=$(getPOPsRemote "${rst_root}")
        local_rst=$(getPOPsLocal "${loc_root}")

    elif [[ "x${sim_name}" == "xtagCH4" ]]; then
        remote_rst=$(getTagCH4Remote "${rst_root}")
        local_rst=$(getTagCH4Local "${loc_root}")

    elif [[ "x${sim_name}" == "xtagCO" ]]; then
        remote_rst=$(getTagCORemote "${rst_root}")
        local_rst=$(getTagCOLocal "${loc_root}")

    elif [[ "x${sim_name}" == "xtagO3" ]]; then
        remote_rst=$(getTagO3Remote "${rst_root}")
        local_rst=$(getTagO3Local "${loc_root}")

    elif [[ "x${sim_name}" == "xTransportTracers" ]]; then
        remote_rst=$(getTracersRemote "${rst_root}")
        local_rst=$(getTracersLocal "${loc_root}")

    fi

    # Copy the restart file to the run directory (for AWS or on a local server)
    copyRemoteToLocal "${remote_rst}" "${local_rst}" "${is_aws}" "${s3_cp}"
    return $?
}
