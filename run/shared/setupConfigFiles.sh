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
#### between GEOS-Chem Classic and GCHP. This script mainly now adds species
#### geoschem_config.yml and modifies diagnostic output based on simulation
#### type.
####
#### Argument: Extra option for full-chemistry simulation (string)
#============================================================================
function set_common_settings() {

    # Check that simulation option is passed
    if [[ $# == 2 ]]; then
        sim_extra_option="${1}"
	model="${2}"
    else
       echo 'Usage: ./setupConfigFiles.sh {sim_extra_option} {model}'
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

	# Change time cycle flag to allow missing species (GCClassic only)
	if [[ "x${model}" == "xGCClassic" ]]; then
	    sed_ie 's|EFYO|CYS|' HEMCO_Config.rc
	fi

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

        # Remove @ from HISTORY diagnostic fields & uncomment default collection
        sed_ie 's|@||'                                 HISTORY.rc
        sed_ie "s|#'Default|'Default|"                 HISTORY.rc

    fi

    #------------------------------------------------------------------------
    # Complex SOA settings
    #------------------------------------------------------------------------
    if [[ "x${sim_extra_option}" == "xbenchmark" ]] || \
       [[ ${sim_extra_option}    =~ "complexSOA" ]] || \
       [[ "x${sim_extra_option}" == "xAPM"       ]]; then

	# Add complex SOA species ASOA* and ASOG* following ALK4
        prev_line='      - ALK4'
        new_line='\      - ASOA1\
      - ASOA2\
      - ASOA3\
      - ASOAN\
      - ASOG1\
      - ASOG2\
      - ASOG3
'
        insert_text "${prev_line}" "${new_line}" geoschem_config.yml

	# Add complex SOA species TSOA* and TSOG* following TOLU
        prev_line='      - TOLU'
        new_line='\      - TSOA0\
      - TSOA1\
      - TSOA2\
      - TSOA3\
      - TSOG0\
      - TSOG1\
      - TSOG2\
      - TSOG3
'
        insert_text "${prev_line}" "${new_line}" geoschem_config.yml

	sed_ie 's/@//' HISTORY.rc
    fi

    # For complexSOA only, remove SOAP and SOAS species from geoschem_config.yml
    if [[ ${sim_extra_option} =~ "complexSOA" ]]; then
	remove_text '      - SOAP' geoschem_config.yml
	remove_text '      - SOAS' geoschem_config.yml
    fi

    #------------------------------------------------------------------------
    # Semivolatile POA settings
    #------------------------------------------------------------------------
    if [[ "x${sim_extra_option}" == "xcomplexSOA_SVPOA" ]]; then

	# Remove non-SVPOA species from geoschem_config.yml
        remove_text '      - OCPI' geoschem_config.yml
        remove_text '      - OCPO' geoschem_config.yml

        # Add semivolatile POA species in geoschem_config.yml
        prev_line='      - N2O5'
        new_line='\      - NAP'
	insert_text "${prev_line}" "${new_line}" geoschem_config.yml

	# Add OPOA* and OPOG* species following OIO
        prev_line='      - OIO'
        new_line='      - OPOA1\
      - OPOA2\
      - OPOG1\
      - OPOG2
'
	insert_text "${prev_line}" "${new_line}" geoschem_config.yml

	# Add POA* and POG* species following PIP
        prev_line='      - PIP'
        new_line='\      - POA1\
      - POA2\
      - POG1\
      - POG2
'
	insert_text "${prev_line}" "${new_line}" geoschem_config.yml

	# Remove the @ from HISTORY.rc diagnostic fields
	sed_ie 's/@//' HISTORY.rc
    fi

    #------------------------------------------------------------------------
    # Acid uptake settings
    #------------------------------------------------------------------------
    if [[ "x${sim_extra_option}" == "xaciduptake" ]]; then

        # Add DSTAL* species after DST4
        prev_line='      - DST4'
        new_line='\      - DSTAL1\
      - DSTAL2\
      - DSTAL3\
      - DSTAL4
'
	insert_text "${prev_line}" "${new_line}" geoschem_config.yml

	# Add NITD* species after NITs.  NOTE: This is non-alphabetical,
	# but avoids double-adding these species after NIT and NITs.
        prev_line='      - NITs'
        new_line='\      - NITD1\
      - NITD2\
      - NITD3\
      - NITD4
'
	insert_text "${prev_line}" "${new_line}" geoschem_config.yml

	# Add SO4* species after SO4s.  NOTE: This is non-alphabetical,
	# but avoids double-adding these species after SO4 and SO4s.
        prev_line='      - SO4s'
        new_line='\      - SO4D1\
      - SO4D2\
      - SO4D3\
      - SO4D4
'
	insert_text "${prev_line}" "${new_line}" geoschem_config.yml

	# Remove the @ from HISTORY.rc diagnostic fields
	sed_ie 's/@//' HISTORY.rc
    fi

    #------------------------------------------------------------------------
    # Marine POA settings
    #------------------------------------------------------------------------
    if [[ "x${sim_extra_option}" == "xmarinePOA" ]]; then

        # Add MOP* species following MONITU
        prev_line='      - MONITU'
        new_line='\      - MOPI\
      - MOPO
'
	insert_text "${prev_line}" "${new_line}" geoschem_config.yml

	# Remove the @ from HISTORY.rc diagnostic fields
	sed_ie 's|@||' HISTORY.rc
    fi

    #------------------------------------------------------------------------
    # RRTMG settings
    #------------------------------------------------------------------------
    if [[ "x${sim_extra_option}" == "xRRTMG" ]]; then

	# Remove @ from HISTORY diagnostic fields & uncomment RRTMG collection
	sed_ie 's|@||'                                 HISTORY.rc
        sed_ie "s|##'RRTMG'|'RRTMG'|"                  HISTORY.rc

	# Issue a warning
	printf "\nWARNING: All RRTMG run options are enabled which will significantly slow down the model!"
        printf "\nEdit geoschem_config.yml and HISTORY.rc in your new run directory to customize options to only"
        printf "\nwhat you need.\n"
    fi

    #------------------------------------------------------------------------
    # TOMAS settings
    #------------------------------------------------------------------------
    if [[ ${sim_extra_option} =~ "TOMAS" ]]; then

	# Remove extra species in extension settings for TOMAS15 simulations
	if [[ "x${sim_extra_option}" == "xTOMAS15" ]]; then
	    sed_ie 's|\/SS16\/SS17\/SS18\/SS19\/SS20\/SS21\/SS22\/SS23\/SS24\/SS25\/SS26\/SS27\/SS28\/SS29\/SS30\/SS31\/SS32\/SS33\/SS34\/SS35\/SS36\/SS37\/SS38\/SS39\/SS40||' HEMCO_Config.rc
	    sed_ie 's|\/DUST16\/DUST17\/DUST18\/DUST19\/DUST20\/DUST21\/DUST22\/DUST23\/DUST24\/DUST25\/DUST26\/DUST27\/DUST28\/DUST29\/DUST30\/DUST31\/DUST32\/DUST33\/DUST34\/DUST35\/DUST36\/DUST37\/DUST38\/DUST39\/DUST40||' HEMCO_Config.rc
	fi

	# Add TOMAS species for the first 15 bins following XYLE
        prev_line='      - XYLE'
        new_line='\      - H2SO4\
      - NK1\
      - NK2\
      - NK3\
      - NK4\
      - NK5\
      - NK6\
      - NK7\
      - NK8\
      - NK9\
      - NK10\
      - NK11\
      - NK12\
      - NK13\
      - NK14\
      - NK15\
      - SF1\
      - SF2\
      - SF3\
      - SF4\
      - SF5\
      - SF6\
      - SF7\
      - SF8\
      - SF9\
      - SF10\
      - SF11\
      - SF12\
      - SF13\
      - SF14\
      - SF15\
      - SS1\
      - SS2\
      - SS3\
      - SS4\
      - SS5\
      - SS6\
      - SS7\
      - SS8\
      - SS9\
      - SS10\
      - SS11\
      - SS12\
      - SS13\
      - SS14\
      - SS15\
      - ECOB1\
      - ECOB2\
      - ECOB3\
      - ECOB4\
      - ECOB5\
      - ECOB6\
      - ECOB7\
      - ECOB8\
      - ECOB9\
      - ECOB10\
      - ECOB11\
      - ECOB12\
      - ECOB13\
      - ECOB14\
      - ECOB15\
      - ECIL1\
      - ECIL2\
      - ECIL3\
      - ECIL4\
      - ECIL5\
      - ECIL6\
      - ECIL7\
      - ECIL8\
      - ECIL9\
      - ECIL10\
      - ECIL11\
      - ECIL12\
      - ECIL13\
      - ECIL14\
      - ECIL15\
      - OCOB1\
      - OCOB2\
      - OCOB3\
      - OCOB4\
      - OCOB5\
      - OCOB6\
      - OCOB7\
      - OCOB8\
      - OCOB9\
      - OCOB10\
      - OCOB11\
      - OCOB12\
      - OCOB13\
      - OCOB14\
      - OCOB15\
      - OCIL1\
      - OCIL2\
      - OCIL3\
      - OCIL4\
      - OCIL5\
      - OCIL6\
      - OCIL7\
      - OCIL8\
      - OCIL9\
      - OCIL10\
      - OCIL11\
      - OCIL12\
      - OCIL13\
      - OCIL14\
      - OCIL15\
      - DUST1\
      - DUST2\
      - DUST3\
      - DUST4\
      - DUST5\
      - DUST6\
      - DUST7\
      - DUST8\
      - DUST9\
      - DUST10\
      - DUST11\
      - DUST12\
      - DUST13\
      - DUST14\
      - DUST15\
      - AW1\
      - AW2\
      - AW3\
      - AW4\
      - AW5\
      - AW6\
      - AW7\
      - AW8\
      - AW9\
      - AW10\
      - AW11\
      - AW12\
      - AW13\
      - AW14\
      - AW15
'
	insert_text "${prev_line}" "${new_line}" geoschem_config.yml

	# Special handling for TOMAS-40 bin simulations
        if [[ ${sim_extra_option} = "TOMAS40" ]]; then

	    # Add NK16-NK40
    	    prev_line='      - NK15'
    	    new_line='\      - NK16\
      - NK17\
      - NK18\
      - NK19\
      - NK20\
      - NK21\
      - NK22\
      - NK23\
      - NK24\
      - NK25\
      - NK26\
      - NK27\
      - NK28\
      - NK29\
      - NK30\
      - NK31\
      - NK32\
      - NK33\
      - NK34\
      - NK35\
      - NK36\
      - NK37\
      - NK38\
      - NK39\
      - NK40
'
	    insert_text "${prev_line}" "${new_line}" geoschem_config.yml

	    # Add SF16-SF40
	    prev_line='      - SF15'
    	    new_line='\      - SF16\
      - SF17\
      - SF18\
      - SF19\
      - SF20\
      - SF21\
      - SF22\
      - SF23\
      - SF24\
      - SF25\
      - SF26\
      - SF27\
      - SF28\
      - SF29\
      - SF30\
      - SF31\
      - SF32\
      - SF33\
      - SF34\
      - SF35\
      - SF36\
      - SF37\
      - SF38\
      - SF39\
      - SF40
'
	    insert_text "${prev_line}" "${new_line}" geoschem_config.yml

	    # Add SS16-SS40
    	    prev_line='      - SS15'
    	    new_line='\      - SS16\
      - SS17\
      - SS18\
      - SS19\
      - SS20\
      - SS21\
      - SS22\
      - SS23\
      - SS24\
      - SS25\
      - SS26\
      - SS27\
      - SS28\
      - SS29\
      - SS30\
      - SS31\
      - SS32\
      - SS33\
      - SS34\
      - SS35\
      - SS36\
      - SS37\
      - SS38\
      - SS39\
      - SS40
'
	    insert_text "${prev_line}" "${new_line}" geoschem_config.yml

	    # Add ECOB16-ECOB40
    	    prev_line='      - ECOB15'
    	    new_line='\      - ECOB16\
      - ECOB17\
      - ECOB18\
      - ECOB19\
      - ECOB20\
      - ECOB21\
      - ECOB22\
      - ECOB23\
      - ECOB24\
      - ECOB25\
      - ECOB26\
      - ECOB27\
      - ECOB28\
      - ECOB29\
      - ECOB30\
      - ECOB31\
      - ECOB32\
      - ECOB33\
      - ECOB34\
      - ECOB35\
      - ECOB36\
      - ECOB37\
      - ECOB38\
      - ECOB39\
      - ECOB40
'
	    insert_text "${prev_line}" "${new_line}" geoschem_config.yml

	    # Add ECIL15-ECIL40
    	    prev_line='      - ECIL15'
    	    new_line='\      - ECIL16\
      - ECIL17\
      - ECIL18\
      - ECIL19\
      - ECIL20\
      - ECIL21\
      - ECIL22\
      - ECIL23\
      - ECIL24\
      - ECIL25\
      - ECIL26\
      - ECIL27\
      - ECIL28\
      - ECIL29\
      - ECIL30\
      - ECIL31\
      - ECIL32\
      - ECIL33\
      - ECIL34\
      - ECIL35\
      - ECIL36\
      - ECIL37\
      - ECIL38\
      - ECIL39\
      - ECIL40
'
	    insert_text "${prev_line}" "${new_line}" geoschem_config.yml

	    # Add OCOB15-OCOB40
    	    prev_line='      - OCOB15'
    	    new_line='\      - OCOB16\
      - OCOB17\
      - OCOB18\
      - OCOB19\
      - OCOB20\
      - OCOB21\
      - OCOB22\
      - OCOB23\
      - OCOB24\
      - OCOB25\
      - OCOB26\
      - OCOB27\
      - OCOB28\
      - OCOB29\
      - OCOB30\
      - OCOB31\
      - OCOB32\
      - OCOB33\
      - OCOB34\
      - OCOB35\
      - OCOB36\
      - OCOB37\
      - OCOB38\
      - OCOB39\
      - OCOB40
'
	    insert_text "${prev_line}" "${new_line}" geoschem_config.yml

	    # Add OCIL16-OCIL40
    	    prev_line='      - OCIL15'
    	    new_line='\      - OCIL16\
      - OCIL17\
      - OCIL18\
      - OCIL19\
      - OCIL20\
      - OCIL21\
      - OCIL22\
      - OCIL23\
      - OCIL24\
      - OCIL25\
      - OCIL26\
      - OCIL27\
      - OCIL28\
      - OCIL29\
      - OCIL30\
      - OCIL31\
      - OCIL32\
      - OCIL33\
      - OCIL34\
      - OCIL35\
      - OCIL36\
      - OCIL37\
      - OCIL38\
      - OCIL39\
      - OCIL40
'
	    insert_text "${prev_line}" "${new_line}" geoschem_config.yml

	    # Add DUST15-DUST40
    	    prev_line='      - DUST15'
    	    new_line='\      - DUST16\
      - DUST17\
      - DUST18\
      - DUST19\
      - DUST20\
      - DUST21\
      - DUST22\
      - DUST23\
      - DUST24\
      - DUST25\
      - DUST26\
      - DUST27\
      - DUST28\
      - DUST29\
      - DUST30\
      - DUST31\
      - DUST32\
      - DUST33\
      - DUST34\
      - DUST35\
      - DUST36\
      - DUST37\
      - DUST38\
      - DUST39\
      - DUST40
'
	    insert_text "${prev_line}" "${new_line}" geoschem_config.yml

	    # Add AW15-AW40
    	    prev_line='      - AW15'
    	    new_line='\      - AW16\
      - AW17\
      - AW18\
      - AW19\
      - AW20\
      - AW21\
      - AW22\
      - AW23\
      - AW24\
      - AW25\
      - AW26\
      - AW27\
      - AW28\
      - AW29\
      - AW30\
      - AW31\
      - AW32\
      - AW33\
      - AW34\
      - AW35\
      - AW36\
      - AW37\
      - AW38\
      - AW39\
      - AW40
'
	    insert_text "${prev_line}" "${new_line}" geoschem_config.yml
        fi

	# Remove the @ from HISTORY.rc diagnostic fields
	sed_ie 's/@//' HISTORY.rc
    fi

    #------------------------------------------------------------------------
    # APM settings
    #------------------------------------------------------------------------
    if [[ "x${sim_extra_option}" == "xAPM" ]]; then

        # Add APM species following XYLE
        prev_line='      - XYLE'
        new_line='\      - APMBCBIN01\
      - APMBCBIN02\
      - APMBCBIN03\
      - APMBCBIN04\
      - APMBCBIN05\
      - APMBCBIN06\
      - APMBCBIN07\
      - APMBCBIN08\
      - APMBCBIN09\
      - APMBCBIN10\
      - APMBCBIN11\
      - APMBCBIN12\
      - APMBCBIN13\
      - APMBCBIN14\
      - APMBCBIN15\
      - APMCTBC1\
      - APMCTBC2\
      - APMCTDST1\
      - APMCTDST2\
      - APMCTOC1\
      - APMCTOC2\
      - APMCTSEA1\
      - APMCTSEA2\
      - APMDSTBIN01\
      - APMDSTBIN02\
      - APMDSTBIN03\
      - APMDSTBIN04\
      - APMDSTBIN05\
      - APMDSTBIN06\
      - APMDSTBIN07\
      - APMDSTBIN08\
      - APMDSTBIN09\
      - APMDSTBIN10\
      - APMDSTBIN11\
      - APMDSTBIN12\
      - APMDSTBIN13\
      - APMDSTBIN14\
      - APMDSTBIN15\
      - APMH2SO4\
      - APMLVSOA\
      - APMLVSOG\
      - APMOCBIN01\
      - APMOCBIN02\
      - APMOCBIN03\
      - APMOCBIN04\
      - APMOCBIN05\
      - APMOCBIN06\
      - APMOCBIN07\
      - APMOCBIN08\
      - APMOCBIN09\
      - APMOCBIN10\
      - APMOCBIN11\
      - APMOCBIN12\
      - APMOCBIN13\
      - APMOCBIN14\
      - APMOCBIN15\
      - APMSEABIN01\
      - APMSEABIN02\
      - APMSEABIN03\
      - APMSEABIN04\
      - APMSEABIN05\
      - APMSEABIN06\
      - APMSEABIN07\
      - APMSEABIN08\
      - APMSEABIN09\
      - APMSEABIN10\
      - APMSEABIN11\
      - APMSEABIN12\
      - APMSEABIN13\
      - APMSEABIN14\
      - APMSEABIN15\
      - APMSEABIN16\
      - APMSEABIN17\
      - APMSEABIN18\
      - APMSEABIN19\
      - APMSEABIN20\
      - APMSPBIN01\
      - APMSPBIN02\
      - APMSPBIN03\
      - APMSPBIN04\
      - APMSPBIN05\
      - APMSPBIN06\
      - APMSPBIN07\
      - APMSPBIN08\
      - APMSPBIN09\
      - APMSPBIN10\
      - APMSPBIN11\
      - APMSPBIN12\
      - APMSPBIN13\
      - APMSPBIN14\
      - APMSPBIN15\
      - APMSPBIN16\
      - APMSPBIN17\
      - APMSPBIN18\
      - APMSPBIN19\
      - APMSPBIN20\
      - APMSPBIN21\
      - APMSPBIN22\
      - APMSPBIN23\
      - APMSPBIN24\
      - APMSPBIN25\
      - APMSPBIN26\
      - APMSPBIN27\
      - APMSPBIN28\
      - APMSPBIN29\
      - APMSPBIN30\
      - APMSPBIN31\
      - APMSPBIN32\
      - APMSPBIN33\
      - APMSPBIN34\
      - APMSPBIN35\
      - APMSPBIN36\
      - APMSPBIN37\
      - APMSPBIN38\
      - APMSPBIN39\
      - APMSPBIN40
'
	insert_text "${prev_line}" "${new_line}" geoschem_config.yml

	# Remove the @ from HISTORY.rc diagnostic fields
	sed_ie 's/@//' HISTORY.rc
    fi
}
#EOC
