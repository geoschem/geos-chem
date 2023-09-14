#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: single_carbon_species.sh
#
# !DESCRIPTION: Updates carbon simulation configuration files so that
#  a simulation with a single species can be performed.
#
# !CALLING SEQUENCE:
#  ./single_carbon
#
# !REVISION HISTORY:
#  06 Jan 2015 - R. Yantosca - Initial version
#  See the subsequent Git history with the gitk browser!
#EOP
#------------------------------------------------------------------------------
#BOC

# List of all advected species in the carbon simulation
ALL_SPECIES=(CO2 CO CH4 OCS)


function isItemInList() {

    #=========================================================================
    # Test if an item is in a list.
    #
    # Arguments:
    # ${1}: The item
    # ${2}: Thie list
    #
    # Returns (via $?)
    # 0 if item is in the list
    # 1 if item is not in the list
    #
    # See stackoverflow.com/questions/8063228/check-if-a-variable-exists-in-a-list-in-bash
    #=========================================================================
    echo "${2}" | tr ' '  '\n' | grep -F -x -q "${1}"
}


function keyValueUpdate() {

    #=========================================================================
    # Runs an sed command to update a "key: value" pair in a file
    #
    # Arguments:
    # ${1} : Key
    # ${2} : Value
    # ${3} : Replacement value
    # ${4) : File in which the text is found
    #=========================================================================
    cmd="s/${1}: ${2}/${1}: ${3}/"
    sed -i -e "$cmd" "${4}"
}


function speciesToExclude() {

    #========================================================================
    # Returns the list of species to exclude
    #
    # Arguments:
    # ${1} : Species to retain
    #
    # Returns:
    # List of species to exclude
    #========================================================================

    # All species
    list=("${ALL_SPECIES[@]}")

    # Keep all species except for the 1st argument
    for i in "${!list[@]}"; do
        if [[ "x${list[i]}" == "x${1}" ]]; then
            unset 'list[i]'
        fi
    done

    # Print the result (that's how we return strings)
    result=""
    for spc in "${list[@]}"; do
        result+="${spc} "
    done
    echo $result
}


function updateGeosChemConfig() {

    #========================================================================
    # Removes advected species from geoschem_config.yml
    #
    # Arguments:
    # ${1} : List of species to exclude
    # ${2} : Path to the run directory
    #========================================================================

    # File to modify
    file="${2}/geoschem_config.yml"

    # Remove advected species (use \< and \> for exact match in order
    # to prevent inadvertently removing CO2 with when spc is CO.
    for spc in ${1}; do
	cmd="/\s*\-\s*\<${spc}\>/d"
        sed -i -e "${cmd}" "${file}"
    done

    # NOTE: CH4 options are already deactivated
    # in the out-of-the-box geoschem_config.yml

#   # If CO is in the exclude list, turn off CO options
#    isItemInList "CO2" "${1}"
#    if [[ $? == 0 ]]; then
#        keys=("use_fullchem_PCO_from_CH4"   \
#              "use_fullchem_PCO_from_NMVOC")
#        for key in ${keys[@]}; do
#            keyValueUpdate "${key}" "true" "false" "${file}"
#        done
#    fi


    # If CO2 is in the exclude list, turn off CO2 options
    isItemInList "CO2" "${1}"
    if [[ $? == 0 ]]; then
        keys=("fossil_fuel_emissions"          \
              "ocean_exchange"                 \
              "balanced_biosphere_exchange"    \
              "net_terrestrial_exchange"       \
              "ship_emissions"                 \
              "aviation_emissions"             \
              "3D_chemical_oxidation_source"   \
              "save_fossil_fuel_in_background" \
              "tag_bio_and_ocean_CO2"          \
              "tag_land_fossil_fuel_CO2"       \
              "tag_global_ship_CO2"            \
              "tag_global_aircraft_CO2"       )
        for key in ${keys[@]}; do
            keyValueUpdate "${key}" "true" "false" "${file}"
        done
    fi
}


function updateHemcoConfig() {

    #========================================================================
    # Removes advected species from geoschem_config.yml
    #
    # Arguments:
    # ${1} : List of species to exclude
    # ${2} : Path to the run directory
    #========================================================================

    # File to be modified
    file="${2}/HEMCO_Config.rc"

    # True/False values in the HEMCO extension switches section
    true="      true"
    false="      false"

    # If CO2 is in the exclude list, turn off CO2 options
    isItemInList "CO2" "${1}"
    if [[ $? == 0 ]]; then
	key="--> USE_CO2_DATA           "
	keyValueUpdate "${key}" "${true}" "${false}" "${file}"
    fi

    # If CO is in the exclude list, turn off CO options
    isItemInList "CO" "${1}"
    if [[ $? == 0 ]]; then
        key="--> USE_CO_DATA            "
	keyValueUpdate "${key}" "${true}" "${false}" "${file}"
    fi

    # If CH4 is in the exclude list, turn off CH4 options
    isItemInList "CH4" "${1}"
    if [[ $? == 0 ]]; then
	key="--> USE_CH4_DATA           "
	keyValueUpdate "${key}" "${true}" "${false}" "${file}"
    fi


    # If OCS is in the exclude list, turn off OCS options
    isItemInList "OCS" "${1}"
    if [[ $? == 0 ]]; then
        key="--> USE_OCS_DATA           "
	keyValueUpdate "${key}" "${true}" "${false}" "${file}"
    fi
}


function updateHemcoDiagn() {

    #========================================================================
    # Comments out lines for unused species in HEMCO_Diagn.rc
    #
    # Arguments:
    # ${1} : List of species to exclude
    # ${2} : Path to the run directory
    #========================================================================

    # File to modify
    file="${2}/HEMCO_Diagn.rc"

    # Remove entries for excluded species
    exclude=("${1}")
    for spc in ${exclude[@]}; do
	sed -i "/Emis${spc}_/d" "${file}"
    done
    sed -i "/#####/d" "${file}"
}


function updateHistory() {

    #========================================================================
    # Removes entries in HISTORY.rc for unused species
    #
    # Arguments:
    # ${1} : List of species to exclude
    # ${2} : Path to the run directory
    #========================================================================

    # File to be modified
    file="${2}/HISTORY.rc"

    # For GCHP: remove entries for species to be excluded
    for spc in ${ALL_SPECIES[@]}; do
	isItemInList "${spc}" "${1}"
	if [[ $? == 0 ]]; then
	    sed -i "/\_${spc}/d"   "${file}"
	    sed -i "/Emis${spc}/d" "${file}"
	fi
    done
}


function updateExtData() {
    
    #========================================================================
    # Removes entries in ExtData.rc for unused species (GCHP only)
    #
    # Arguments:
    # ${1} : List of species to exclude
    # ${2} : Path to the run directory
    #========================================================================

    # Skip if there is no ExtData.rc file
    [[ ! -f "${2}/ExtData.rc" ]] && return 0
    
}


function singleCarbonSpecies() {

    #========================================================================
    # Main function
    #
    # Arguments:
    # ${1} : Species that you wish to retain
    # ${2} : Path to the run directory    
    #========================================================================

    # Error check arguments
    if [[ "x${1}" == "x" ]]; then
        echo "Need to pass the species to retain!"
        exit 1
    fi

    # Path to the run directory
    if [[ "x${2}" == "x" ]]; then
	rundir="."
    else
	rundir="${2}"
    fi
    
    # Get species to exclude
    exclude=$(speciesToExclude "${1}")
    
    # Update configuration files
    updateGeosChemConfig "${exclude}" "${rundir}"
    updateHemcoConfig    "${exclude}" "${rundir}"
    updateHemcoDiagn     "${exclude}" "${rundir}"
    updateHistory        "${exclude}" "${rundir}"
    updateExtData        "${exclude}" "${rundir}"
}
