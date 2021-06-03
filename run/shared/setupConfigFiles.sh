#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: setupConfigFiles.sh
#
# !DESCRIPTION: Defines utility functions to update configuration file
#  settings when creating a run directory.
#\\
#\\
# !REVISION HISTORY:
#  Initial version: M. Sulprizio, 6/24/2020
#  See the subsequent Git history with the gitk browser!
#EOP
#------------------------------------------------------------------------------
#BOC

#=============================================================================
#### Replacement for `sed -i -e` that works on both MacOS and Linux
#=============================================================================
function sed_ie() {
    REGEX=${1}
    FILE=${2}
    if [[ "x$(uname -s)" == "xDarwin" ]]; then
	sed -i '' -e "${REGEX}" "${FILE}"          # MacOS/Darwin
    else
	sed -i -e "${REGEX}" "${FILE}"             # GNU/Linux
    fi
}

#=============================================================================
#### Define function to replace values in config files
#=============================================================================
function replace_colon_sep_val() {
    KEY=${1}
    VALUE=${2}
    FILE=${3}

    # Debug print (leave commented out)
    # printf '%-30s : %-20s %-20s\n' "${KEY//\\}" "${VALUE}" "${FILE}"

    # Replace value in line starting with 'whitespace + key + whitespace + : +
    # whitespace + value' where whitespace is variable length including none
    #
    # MacOS sed does not allow you to use \t for tab.  The quick fix is
    # to use printf to save a tab character to a variable, and to use
    # that in the regular expression everywhere you would have used \t.
    # See the Github issue geoschem/geos-chem #617. (bmy, 2/23/21)
    TAB=$(printf "\t")
    REGEX="s|^\([${TAB} ]*${KEY}[${TAB} ]*:[${TAB} ]*\).*|\1${VALUE}|"
    sed_ie "${REGEX}" "${FILE}"
}

#============================================================================
#### Define function to remove line(s) in config files
#============================================================================
function remove_text() {
    REGEX="/${1}/d"
    FILE=${2}
    sed_ie "${REGEX}" "${FILE}"
}

#============================================================================
#### Define function to insert text between two lines of a config file
#============================================================================
function insert_text() {
    PREV_LINE=${1}
    NEW_LINE=${2}
    FILE=${3}

    # The POSIX standard (which works on both MacOS & GNU/Linux) calls
    # for splitting the regular expression with a hard return instead
    # of using a "\n" newline character.
    REGEX="/${PREV_LINE}/a \\
${NEW_LINE}"
    sed_ie "${REGEX}" "${FILE}"
}

#============================================================================
#### Define function to update config file default settings based on
#### simulation selected. All settings changed in this function are common
#### between GEOS-Chem Classic and GCHP.
####
#### Argument: Extra option for full-chemistry simulation (string)
#============================================================================
function set_common_settings() {

    # Check that simulation option is passed
    if [[ $# == 1 ]]; then
        sim_extra_option=${1}
    else
       echo 'Usage: ./setupConfigFiles.sh {sim_extra_option}'
       exit 1
    fi

    valid_options=( 'standard' 'benchmark' 'complexSOA' 'complexSOA_SVPOA' \
                    'aciduptake' 'marinePOA' 'TOMAS15' 'TOMAS40' 'APM' 'RRTMG' )
    for i in "${valid_options[@]}"; do
        if [ "$i" == "$yourValue" ] ; then
            echo "Found"
        fi
    done

    #------------------------------------------------------------------------
    # Benchmark settings
    #------------------------------------------------------------------------
    if [[ "x${sim_extra_option}" == "xbenchmark" ]]; then
        replace_colon_sep_val '--> OFFLINE_DUST'        false HEMCO_Config.rc
        replace_colon_sep_val '--> OFFLINE_BIOGENICVOC' false HEMCO_Config.rc
        replace_colon_sep_val '--> OFFLINE_SEASALT'     false HEMCO_Config.rc
        replace_colon_sep_val '--> OFFLINE_SOILNOX'     false HEMCO_Config.rc
        sed_ie 's|DustDead               : off|DustDead               : on |' HEMCO_Config.rc
        sed_ie 's|SoilNOx                : off|SoilNOx                : on |' HEMCO_Config.rc
        sed_ie 's|SeaSalt                : off|SeaSalt                : on |' HEMCO_Config.rc
        sed_ie 's|NO     0      3 |NO     104    -1|' HEMCO_Diagn.rc   # Use online soil NOx (ExtNr=104)
	sed_ie 's|SALA  0      3 |SALA  107    -1|'   HEMCO_Diagn.rc   # Use online sea salt (ExtNr=107)
	sed_ie 's|SALC  0      3 |SALC  107    -1|'   HEMCO_Diagn.rc   #   "   "
	sed_ie 's|AL  0      3 |AL  107    -1|'       HEMCO_Diagn.rc   #   "   "
	sed_ie 's|CL  0      3 |CL  107    -1|'       HEMCO_Diagn.rc   #   "   "
        sed_ie 's|0      3 |105    -1|'               HEMCO_Diagn.rc   # Use online dust (ExtNr=105)
        sed_ie 's|0      4 |108    -1|'               HEMCO_Diagn.rc   # Use MEGAN (ExtNr=108)
        sed_ie 's|NH3    105    -1|NH3    0      3 |' HEMCO_Diagn.rc   # NaturalNH3 is always ExtNr=0
        sed_ie 's|ALD2   105    -1|ALD2   0      3 |' HEMCO_Diagn.rc   # PlantDecay is always ExtNr=0
        sed_ie 's|EOH    105    -1|EOH    0      3 |' HEMCO_Diagn.rc   # PlantDecay is always ExtNr=0
        sed_ie 's|#Inv|Inv|'                          HEMCO_Diagn.rc

	# Turn @ into # characters for the benchmark simulation,
	# which should cause MAPL to skip reading these lines.
	# This is a workaround for a "input file to long" MAPL error.
	sed_ie 's|@|#|'                                HISTORY.rc

	# Remove the first comment character on diagnostics
        sed_ie "s|#'|'|"                               HISTORY.rc
    fi

    #------------------------------------------------------------------------
    # Standard settings
    #------------------------------------------------------------------------
    if [[ "x${sim_extra_option}" == "xnone" ]]; then
	sed_ie 's/@//' HISTORY.rc
    fi

    #------------------------------------------------------------------------
    # Complex SOA settings
    #------------------------------------------------------------------------
    if [[ "x${sim_extra_option}" == "xbenchmark" ]] || \
       [[ ${sim_extra_option}    =~ "complexSOA" ]] || \
       [[ "x${sim_extra_option}" == "xAPM"       ]]; then

        # Turn on complex SOA option in input.geos
        replace_colon_sep_val 'Online COMPLEX SOA' T input.geos

	# Add complex SOA species ASOA* and ASOG* following ALK4
        prev_line='Species name            : ALK4'
        new_line='\Species name            : ASOA1\
Species name            : ASOA2\
Species name            : ASOA3\
Species name            : ASOAN\
Species name            : ASOG1\
Species name            : ASOG2\
Species name            : ASOG3
'
        insert_text "${prev_line}" "${new_line}" input.geos

	# Add complex SOA species TSOA* and TSOG* following TOLU
        prev_line='Species name            : TOLU'
        new_line='\Species name            : TSOA0\
Species name            : TSOA1\
Species name            : TSOA2\
Species name            : TSOA3\
Species name            : TSOG0\
Species name            : TSOG1\
Species name            : TSOG2\
Species name            : TSOG3
'
        insert_text "${prev_line}" "${new_line}" input.geos

	sed_ie 's/@//' HISTORY.rc
    fi

    # For complexSOA only, remove SOAP and SOAS species from input.geos
    if [[ ${sim_extra_option} =~ "complexSOA" ]]; then
	remove_text 'Species name            : SOAP' input.geos
	remove_text 'Species name            : SOAS' input.geos
    fi

    #------------------------------------------------------------------------
    # Semivolatile POA settings
    #------------------------------------------------------------------------
    if [[ "x${sim_extra_option}" == "xcomplexSOA_SVPOA" ]]; then

        # Turn on semivolatile POA option in input.geos
        replace_colon_sep_val '=> Semivolatile POA?' T input.geos

	# Remove non-SVPOA species from input.geos
        remove_text 'Species name            : OCPI' input.geos
        remove_text 'Species name            : OCPO' input.geos

        # Add semivolatile POA species in input.geos
        prev_line='Species name            : N2O5'
        new_line='\Species name            : NAP'
	insert_text "${prev_line}" "${new_line}" input.geos

	# Add OPOA* and OPOG* species following OIO
        prev_line='Species name            : OIO'
        new_line='Species name            : OPOA1\
Species name            : OPOA2\
Species name            : OPOG1\
Species name            : OPOG2
'
	insert_text "${prev_line}" "${new_line}" input.geos

	# Add POA* and POG* species following PIP
        prev_line='Species name            : PIP'
        new_line='\Species name            : POA1\
Species name            : POA2\
Species name            : POG1\
Species name            : POG2
'
	insert_text "${prev_line}" "${new_line}" input.geos

	# Remove the @ from HISTORY.rc diagnostic fields
	sed_ie 's/@//' HISTORY.rc
    fi

    #------------------------------------------------------------------------
    # Acid uptake settings
    #------------------------------------------------------------------------
    if [[ "x${sim_extra_option}" == "xaciduptake" ]]; then
        replace_colon_sep_val 'DustAlk' on HEMCO_Config.rc
        replace_colon_sep_val '=> Acidic uptake ?' T input.geos

        # Add DSTAL* species after DST4
        prev_line='Species name            : DST4'
        new_line='\Species name            : DSTAL1\
Species name            : DSTAL2\
Species name            : DSTAL3\
Species name            : DSTAL4
'
	insert_text "${prev_line}" "${new_line}" input.geos

	# Add NITD* species after NITs.  NOTE: This is non-alphabetical,
	# but avoids double-adding these species after NIT and NITs.
        prev_line='Species name            : NITs'
        new_line='\Species name            : NITD1\
Species name            : NITD2\
Species name            : NITD3\
Species name            : NITD4
'
	insert_text "${prev_line}" "${new_line}" input.geos

	# Add SO4* species after SO4s.  NOTE: This is non-alphabetical,
	# but avoids double-adding these species after SO4 and SO4s.
        prev_line='Species name            : SO4s'
        new_line='\Species name            : SO4D1\
Species name            : SO4D2\
Species name            : SO4D3\
Species name            : SO4D4
'
	insert_text "${prev_line}" "${new_line}" input.geos

	# Remove the @ from HISTORY.rc diagnostic fields
	sed_ie 's/@//' HISTORY.rc
    fi

    #------------------------------------------------------------------------
    # Marine POA settings
    #------------------------------------------------------------------------
    if [[ "x${sim_extra_option}" == "xmarinePOA" ]]; then
        replace_colon_sep_val 'SeaSalt'                 on HEMCO_Config.rc
        replace_colon_sep_val ' => MARINE ORG AEROSOLS' T  input.geos

        # Add MOP* species following MONITU
        prev_line='Species name            : MONITU'
        new_line='\Species name            : MOPI\
Species name            : MOPO
'
	insert_text "${prev_line}" "${new_line}" input.geos

	# Remove the @ from HISTORY.rc diagnostic fields
	sed_ie 's|@||' HISTORY.rc
    fi

    #------------------------------------------------------------------------
    # RRTMG settings
    #------------------------------------------------------------------------
    if [[ "x${sim_extra_option}" == "xRRTMG" ]]; then

        replace_colon_sep_val 'Turn on RRTMG?'       T input.geos
        replace_colon_sep_val 'Calculate LW fluxes?' T input.geos
        replace_colon_sep_val 'Calculate SW fluxes?' T input.geos
        replace_colon_sep_val 'Clear-sky flux?'      T input.geos
        replace_colon_sep_val 'All-sky flux?'        T input.geos
        replace_colon_sep_val '--> RRTMG'         true HEMCO_Config.rc

	# Remove @ from HISTORY diagnostic fields & uncomment RRTMG collection
	sed_ie 's|@||'                                 HISTORY.rc
        sed_ie "s|##'RRTMG'|'RRTMG'|"                  HISTORY.rc

	# Issue a warning
	printf "\nWARNING: All RRTMG run options are enabled which will significantly slow down the model!"
        printf "\nEdit input.geos and HISTORY.rc in your new run directory to customize options to only"
        printf "\nwhat you need.\n"
    fi

    #------------------------------------------------------------------------
    # TOMAS settings
    #------------------------------------------------------------------------
    if [[ ${sim_extra_option} =~ "TOMAS" ]]; then
        replace_colon_sep_val ' => Use non-local PBL?'   F    input.geos
        replace_colon_sep_val 'Use linear. strat. chem?' F    input.geos
        replace_colon_sep_val '=> Online O3 from model'  F    input.geos

	# Turn on TOMAS size-resolved dust & seasalt extensions
	sed_ie 's/DustDead               : on /DustDead               : off/' HEMCO_Config.rc
	sed_ie 's/TOMAS_DustDead         : off/TOMAS_DustDead         : on /' HEMCO_Config.rc
	sed_ie 's/TOMAS_Jeagle           : off/TOMAS_Jeagle           : on /' HEMCO_Config.rc

	# Remove extra species in extension settings for TOMAS15 simulations
	if [[ "x${sim_extra_option}" == "xTOMAS15" ]]; then
	    sed_ie 's|\/SS16\/SS17\/SS18\/SS19\/SS20\/SS21\/SS22\/SS23\/SS24\/SS25\/SS26\/SS27\/SS28\/SS29\/SS30\/SS31\/SS32\/SS33\/SS34\/SS35\/SS36\/SS37\/SS38\/SS39\/SS40||' HEMCO_Config.rc
	    sed_ie 's|\/DUST16\/DUST17\/DUST18\/DUST19\/DUST20\/DUST21\/DUST22\/DUST23\/DUST24\/DUST25\/DUST26\/DUST27\/DUST28\/DUST29\/DUST30\/DUST31\/DUST32\/DUST33\/DUST34\/DUST35\/DUST36\/DUST37\/DUST38\/DUST39\/DUST40||' HEMCO_Config.rc
	fi

	# Add TOMAS species for the first 15 bins following XYLE
        prev_line='Species name            : XYLE'
        new_line='\Species name            : H2SO4\
Species name            : NK1\
Species name            : NK2\
Species name            : NK3\
Species name            : NK4\
Species name            : NK5\
Species name            : NK6\
Species name            : NK7\
Species name            : NK8\
Species name            : NK9\
Species name            : NK10\
Species name            : NK11\
Species name            : NK12\
Species name            : NK13\
Species name            : NK14\
Species name            : NK15\
Species name            : SF1\
Species name            : SF2\
Species name            : SF3\
Species name            : SF4\
Species name            : SF5\
Species name            : SF6\
Species name            : SF7\
Species name            : SF8\
Species name            : SF9\
Species name            : SF10\
Species name            : SF11\
Species name            : SF12\
Species name            : SF13\
Species name            : SF14\
Species name            : SF15\
Species name            : SS1\
Species name            : SS2\
Species name            : SS3\
Species name            : SS4\
Species name            : SS5\
Species name            : SS6\
Species name            : SS7\
Species name            : SS8\
Species name            : SS9\
Species name            : SS10\
Species name            : SS11\
Species name            : SS12\
Species name            : SS13\
Species name            : SS14\
Species name            : SS15\
Species name            : ECOB1\
Species name            : ECOB2\
Species name            : ECOB3\
Species name            : ECOB4\
Species name            : ECOB5\
Species name            : ECOB6\
Species name            : ECOB7\
Species name            : ECOB8\
Species name            : ECOB9\
Species name            : ECOB10\
Species name            : ECOB11\
Species name            : ECOB12\
Species name            : ECOB13\
Species name            : ECOB14\
Species name            : ECOB15\
Species name            : ECIL1\
Species name            : ECIL2\
Species name            : ECIL3\
Species name            : ECIL4\
Species name            : ECIL5\
Species name            : ECIL6\
Species name            : ECIL7\
Species name            : ECIL8\
Species name            : ECIL9\
Species name            : ECIL10\
Species name            : ECIL11\
Species name            : ECIL12\
Species name            : ECIL13\
Species name            : ECIL14\
Species name            : ECIL15\
Species name            : OCOB1\
Species name            : OCOB2\
Species name            : OCOB3\
Species name            : OCOB4\
Species name            : OCOB5\
Species name            : OCOB6\
Species name            : OCOB7\
Species name            : OCOB8\
Species name            : OCOB9\
Species name            : OCOB10\
Species name            : OCOB11\
Species name            : OCOB12\
Species name            : OCOB13\
Species name            : OCOB14\
Species name            : OCOB15\
Species name            : OCIL1\
Species name            : OCIL2\
Species name            : OCIL3\
Species name            : OCIL4\
Species name            : OCIL5\
Species name            : OCIL6\
Species name            : OCIL7\
Species name            : OCIL8\
Species name            : OCIL9\
Species name            : OCIL10\
Species name            : OCIL11\
Species name            : OCIL12\
Species name            : OCIL13\
Species name            : OCIL14\
Species name            : OCIL15\
Species name            : DUST1\
Species name            : DUST2\
Species name            : DUST3\
Species name            : DUST4\
Species name            : DUST5\
Species name            : DUST6\
Species name            : DUST7\
Species name            : DUST8\
Species name            : DUST9\
Species name            : DUST10\
Species name            : DUST11\
Species name            : DUST12\
Species name            : DUST13\
Species name            : DUST14\
Species name            : DUST15\
Species name            : AW1\
Species name            : AW2\
Species name            : AW3\
Species name            : AW4\
Species name            : AW5\
Species name            : AW6\
Species name            : AW7\
Species name            : AW8\
Species name            : AW9\
Species name            : AW10\
Species name            : AW11\
Species name            : AW12\
Species name            : AW13\
Species name            : AW14\
Species name            : AW15
'
	insert_text "${prev_line}" "${new_line}" input.geos

	# Special handling for TOMAS-40 bin simulations
        if [[ ${sim_extra_option} = "TOMAS40" ]]; then

	    # Add NK16-NK40
    	    prev_line='Species name            : NK15'
    	    new_line='\Species name            : NK16\
Species name            : NK17\
Species name            : NK18\
Species name            : NK19\
Species name            : NK20\
Species name            : NK21\
Species name            : NK22\
Species name            : NK23\
Species name            : NK24\
Species name            : NK25\
Species name            : NK26\
Species name            : NK27\
Species name            : NK28\
Species name            : NK29\
Species name            : NK30\
Species name            : NK31\
Species name            : NK32\
Species name            : NK33\
Species name            : NK34\
Species name            : NK35\
Species name            : NK36\
Species name            : NK37\
Species name            : NK38\
Species name            : NK39\
Species name            : NK40
'
	    insert_text "${prev_line}" "${new_line}" input.geos

	    # Add SF16-SF40
	    prev_line='Species name            : SF15'
    	    new_line='\Species name            : SF16\
Species name            : SF17\
Species name            : SF18\
Species name            : SF19\
Species name            : SF20\
Species name            : SF21\
Species name            : SF22\
Species name            : SF23\
Species name            : SF24\
Species name            : SF25\
Species name            : SF26\
Species name            : SF27\
Species name            : SF28\
Species name            : SF29\
Species name            : SF30\
Species name            : SF31\
Species name            : SF32\
Species name            : SF33\
Species name            : SF34\
Species name            : SF35\
Species name            : SF36\
Species name            : SF37\
Species name            : SF38\
Species name            : SF39\
Species name            : SF40
'
	    insert_text "${prev_line}" "${new_line}" input.geos

	    # Add SS16-SS40
    	    prev_line='Species name            : SS15'
    	    new_line='\Species name            : SS16\
Species name            : SS17\
Species name            : SS18\
Species name            : SS19\
Species name            : SS20\
Species name            : SS21\
Species name            : SS22\
Species name            : SS23\
Species name            : SS24\
Species name            : SS25\
Species name            : SS26\
Species name            : SS27\
Species name            : SS28\
Species name            : SS29\
Species name            : SS30\
Species name            : SS31\
Species name            : SS32\
Species name            : SS33\
Species name            : SS34\
Species name            : SS35\
Species name            : SS36\
Species name            : SS37\
Species name            : SS38\
Species name            : SS39\
Species name            : SS40
'
	    insert_text "${prev_line}" "${new_line}" input.geos

	    # Add ECOB16-ECOB40
    	    prev_line='Species name            : ECOB15'
    	    new_line='\Species name            : ECOB16\
Species name            : ECOB17\
Species name            : ECOB18\
Species name            : ECOB19\
Species name            : ECOB20\
Species name            : ECOB21\
Species name            : ECOB22\
Species name            : ECOB23\
Species name            : ECOB24\
Species name            : ECOB25\
Species name            : ECOB26\
Species name            : ECOB27\
Species name            : ECOB28\
Species name            : ECOB29\
Species name            : ECOB30\
Species name            : ECOB31\
Species name            : ECOB32\
Species name            : ECOB33\
Species name            : ECOB34\
Species name            : ECOB35\
Species name            : ECOB36\
Species name            : ECOB37\
Species name            : ECOB38\
Species name            : ECOB39\
Species name            : ECOB40
'
	    insert_text "${prev_line}" "${new_line}" input.geos

	    # Add ECIL15-ECIL40
    	    prev_line='Species name            : ECIL15'
    	    new_line='\Species name            : ECIL16\
Species name            : ECIL17\
Species name            : ECIL18\
Species name            : ECIL19\
Species name            : ECIL20\
Species name            : ECIL21\
Species name            : ECIL22\
Species name            : ECIL23\
Species name            : ECIL24\
Species name            : ECIL25\
Species name            : ECIL26\
Species name            : ECIL27\
Species name            : ECIL28\
Species name            : ECIL29\
Species name            : ECIL30\
Species name            : ECIL31\
Species name            : ECIL32\
Species name            : ECIL33\
Species name            : ECIL34\
Species name            : ECIL35\
Species name            : ECIL36\
Species name            : ECIL37\
Species name            : ECIL38\
Species name            : ECIL39\
Species name            : ECIL40
'
	    insert_text "${prev_line}" "${new_line}" input.geos

	    # Add OCOB15-OCOB40
    	    prev_line='Species name            : OCOB15'
    	    new_line='\Species name            : OCOB16\
Species name            : OCOB17\
Species name            : OCOB18\
Species name            : OCOB19\
Species name            : OCOB20\
Species name            : OCOB21\
Species name            : OCOB22\
Species name            : OCOB23\
Species name            : OCOB24\
Species name            : OCOB25\
Species name            : OCOB26\
Species name            : OCOB27\
Species name            : OCOB28\
Species name            : OCOB29\
Species name            : OCOB30\
Species name            : OCOB31\
Species name            : OCOB32\
Species name            : OCOB33\
Species name            : OCOB34\
Species name            : OCOB35\
Species name            : OCOB36\
Species name            : OCOB37\
Species name            : OCOB38\
Species name            : OCOB39\
Species name            : OCOB40
'
	    insert_text "${prev_line}" "${new_line}" input.geos

	    # Add OCIL16-OCIL40
    	    prev_line='Species name            : OCIL15'
    	    new_line='\Species name            : OCIL16\
Species name            : OCIL17\
Species name            : OCIL18\
Species name            : OCIL19\
Species name            : OCIL20\
Species name            : OCIL21\
Species name            : OCIL22\
Species name            : OCIL23\
Species name            : OCIL24\
Species name            : OCIL25\
Species name            : OCIL26\
Species name            : OCIL27\
Species name            : OCIL28\
Species name            : OCIL29\
Species name            : OCIL30\
Species name            : OCIL31\
Species name            : OCIL32\
Species name            : OCIL33\
Species name            : OCIL34\
Species name            : OCIL35\
Species name            : OCIL36\
Species name            : OCIL37\
Species name            : OCIL38\
Species name            : OCIL39\
Species name            : OCIL40
'
	    insert_text "${prev_line}" "${new_line}" input.geos

	    # Add DUST15-DUST40
    	    prev_line='Species name            : DUST15'
    	    new_line='\Species name            : DUST16\
Species name            : DUST17\
Species name            : DUST18\
Species name            : DUST19\
Species name            : DUST20\
Species name            : DUST21\
Species name            : DUST22\
Species name            : DUST23\
Species name            : DUST24\
Species name            : DUST25\
Species name            : DUST26\
Species name            : DUST27\
Species name            : DUST28\
Species name            : DUST29\
Species name            : DUST30\
Species name            : DUST31\
Species name            : DUST32\
Species name            : DUST33\
Species name            : DUST34\
Species name            : DUST35\
Species name            : DUST36\
Species name            : DUST37\
Species name            : DUST38\
Species name            : DUST39\
Species name            : DUST40
'
	    insert_text "${prev_line}" "${new_line}" input.geos

	    # Add AW15-AW40
    	    prev_line='Species name            : AW15'
    	    new_line='\Species name            : AW16\
Species name            : AW17\
Species name            : AW18\
Species name            : AW19\
Species name            : AW20\
Species name            : AW21\
Species name            : AW22\
Species name            : AW23\
Species name            : AW24\
Species name            : AW25\
Species name            : AW26\
Species name            : AW27\
Species name            : AW28\
Species name            : AW29\
Species name            : AW30\
Species name            : AW31\
Species name            : AW32\
Species name            : AW33\
Species name            : AW34\
Species name            : AW35\
Species name            : AW36\
Species name            : AW37\
Species name            : AW38\
Species name            : AW39\
Species name            : AW40
'
	    insert_text "${prev_line}" "${new_line}" input.geos
        fi

	# Remove the @ from HISTORY.rc diagnostic fields
	sed_ie 's/@//' HISTORY.rc
    fi

    #------------------------------------------------------------------------
    # APM settings
    #------------------------------------------------------------------------
    if [[ "x${sim_extra_option}" == "xAPM" ]]; then

        # Add APM species following XYLE
        prev_line='Species name            : XYLE'
        new_line='\Species name            : APMBCBIN01\
Species name            : APMBCBIN02\
Species name            : APMBCBIN03\
Species name            : APMBCBIN04\
Species name            : APMBCBIN05\
Species name            : APMBCBIN06\
Species name            : APMBCBIN07\
Species name            : APMBCBIN08\
Species name            : APMBCBIN09\
Species name            : APMBCBIN10\
Species name            : APMBCBIN11\
Species name            : APMBCBIN12\
Species name            : APMBCBIN13\
Species name            : APMBCBIN14\
Species name            : APMBCBIN15\
Species name            : APMCTBC1\
Species name            : APMCTBC2\
Species name            : APMCTDST1\
Species name            : APMCTDST2\
Species name            : APMCTOC1\
Species name            : APMCTOC2\
Species name            : APMCTSEA1\
Species name            : APMCTSEA2\
Species name            : APMDSTBIN01\
Species name            : APMDSTBIN02\
Species name            : APMDSTBIN03\
Species name            : APMDSTBIN04\
Species name            : APMDSTBIN05\
Species name            : APMDSTBIN06\
Species name            : APMDSTBIN07\
Species name            : APMDSTBIN08\
Species name            : APMDSTBIN09\
Species name            : APMDSTBIN10\
Species name            : APMDSTBIN11\
Species name            : APMDSTBIN12\
Species name            : APMDSTBIN13\
Species name            : APMDSTBIN14\
Species name            : APMDSTBIN15\
Species name            : APMH2SO4\
Species name            : APMLVSOA\
Species name            : APMLVSOG\
Species name            : APMOCBIN01\
Species name            : APMOCBIN02\
Species name            : APMOCBIN03\
Species name            : APMOCBIN04\
Species name            : APMOCBIN05\
Species name            : APMOCBIN06\
Species name            : APMOCBIN07\
Species name            : APMOCBIN08\
Species name            : APMOCBIN09\
Species name            : APMOCBIN10\
Species name            : APMOCBIN11\
Species name            : APMOCBIN12\
Species name            : APMOCBIN13\
Species name            : APMOCBIN14\
Species name            : APMOCBIN15\
Species name            : APMSEABIN01\
Species name            : APMSEABIN02\
Species name            : APMSEABIN03\
Species name            : APMSEABIN04\
Species name            : APMSEABIN05\
Species name            : APMSEABIN06\
Species name            : APMSEABIN07\
Species name            : APMSEABIN08\
Species name            : APMSEABIN09\
Species name            : APMSEABIN10\
Species name            : APMSEABIN11\
Species name            : APMSEABIN12\
Species name            : APMSEABIN13\
Species name            : APMSEABIN14\
Species name            : APMSEABIN15\
Species name            : APMSEABIN16\
Species name            : APMSEABIN17\
Species name            : APMSEABIN18\
Species name            : APMSEABIN19\
Species name            : APMSEABIN20\
Species name            : APMSPBIN01\
Species name            : APMSPBIN02\
Species name            : APMSPBIN03\
Species name            : APMSPBIN04\
Species name            : APMSPBIN05\
Species name            : APMSPBIN06\
Species name            : APMSPBIN07\
Species name            : APMSPBIN08\
Species name            : APMSPBIN09\
Species name            : APMSPBIN10\
Species name            : APMSPBIN11\
Species name            : APMSPBIN12\
Species name            : APMSPBIN13\
Species name            : APMSPBIN14\
Species name            : APMSPBIN15\
Species name            : APMSPBIN16\
Species name            : APMSPBIN17\
Species name            : APMSPBIN18\
Species name            : APMSPBIN19\
Species name            : APMSPBIN20\
Species name            : APMSPBIN21\
Species name            : APMSPBIN22\
Species name            : APMSPBIN23\
Species name            : APMSPBIN24\
Species name            : APMSPBIN25\
Species name            : APMSPBIN26\
Species name            : APMSPBIN27\
Species name            : APMSPBIN28\
Species name            : APMSPBIN29\
Species name            : APMSPBIN30\
Species name            : APMSPBIN31\
Species name            : APMSPBIN32\
Species name            : APMSPBIN33\
Species name            : APMSPBIN34\
Species name            : APMSPBIN35\
Species name            : APMSPBIN36\
Species name            : APMSPBIN37\
Species name            : APMSPBIN38\
Species name            : APMSPBIN39\
Species name            : APMSPBIN40
'
	insert_text "${prev_line}" "${new_line}" input.geos

	# Remove the @ from HISTORY.rc diagnostic fields
	sed_ie 's/@//' HISTORY.rc
    fi
}
#EOC
