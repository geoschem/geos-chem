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
#### Define function to replace values in config files
#=============================================================================
replace_colon_sep_val() {
    KEY=$1
    VALUE=$2
    FILE=$3
    #printf '%-30s : %-20s %-20s\n' "${KEY//\\}" "${VALUE}" "${FILE}"

    # replace value in line starting with 'whitespace + key + whitespace + : +
    # whitespace + value' where whitespace is variable length including none
    sed -i -e "s|^\([\t ]*${KEY}[\t ]*:[\t ]*\).*|\1${VALUE}|" ${FILE}
}

#============================================================================
#### Define function to remove line(s) in config files
#============================================================================
remove_text() {
    VALUE=$1
    FILE=$2
    sed -i -e "/${VALUE}/d" ${FILE}
}

#============================================================================
#### Define function to update config file default settings based on
#### simulation selected. All settings changed in this function are common
#### between GEOS-Chem Classic and GCHP.
####
#### Argument: Extra option for full-chemistry simulation (string)
#============================================================================
set_common_settings() {

    # Check that simulation option is passed
    if [[ $# == 1 ]]; then
        sim_extra_option=$1
    else
       echo "Usage: ./setupConfigFiles.sh {sim_extra_option}"
       exit
    fi

    valid_options=( "standard" "benchmark" "complexSOA" "complexSOA_SVPOA" \
                    "aciduptake" "marinePOA" "TOMAS15" "TOMAS40" "APM" "RRTMG" )
    for i in "${valid_options[@]}"; do
        if [ "$i" == "$yourValue" ] ; then
            echo "Found"
        fi
    done

    #------------------------------------------------------------------------
    # Benchmark settings
    #------------------------------------------------------------------------
    if [[ "x${sim_extra_option}" == "xbenchmark" ]]; then
        replace_colon_sep_val "--> OFFLINE_DUST"        false HEMCO_Config.rc
        replace_colon_sep_val "--> OFFLINE_BIOGENICVOC" false HEMCO_Config.rc
        replace_colon_sep_val "--> OFFLINE_SEASALT"     false HEMCO_Config.rc
        replace_colon_sep_val "--> OFFLINE_SOILNOX"     false HEMCO_Config.rc
        sed -i -e "s|DustDead               : off|DustDead               : on |" HEMCO_Config.rc
        sed -i -e "s|SoilNOx                : off|SoilNOx                : on |" HEMCO_Config.rc
        sed -i -e "s|SeaSalt                : off|SeaSalt                : on |" HEMCO_Config.rc
        sed -i -e "s|NO     0      3 |NO     104    -1|" HEMCO_Diagn.rc   # Use online soil NOx (ExtNr=104)
	sed -i -e "s|SALA  0      3 |SALA  107    -1|"   HEMCO_Diagn.rc   # Use online sea salt (ExtNr=107)
	sed -i -e "s|SALC  0      3 |SALC  107    -1|"   HEMCO_Diagn.rc   #   "   "
	sed -i -e "s|AL  0      3 |AL  107    -1|"       HEMCO_Diagn.rc   #   "   "
	sed -i -e "s|CL  0      3 |CL  107    -1|"       HEMCO_Diagn.rc   #   "   "
        sed -i -e "s|0      3 |105    -1|"               HEMCO_Diagn.rc   # Use online dust (ExtNr=105)
        sed -i -e "s|0      4 |108    -1|"               HEMCO_Diagn.rc   # Use MEGAN (ExtNr=108) 
        sed -i -e "s|NH3    105    -1|NH3    0      3 |" HEMCO_Diagn.rc   # NaturalNH3 is always ExtNr=0
        sed -i -e "s|ALD2   105    -1|ALD2   0      3 |" HEMCO_Diagn.rc   # PlantDecay is always ExtNr=0
        sed -i -e "s|EOH    105    -1|EOH    0      3 |" HEMCO_Diagn.rc   # PlantDecay is always ExtNr=0
        sed -i -e "s|#Inv|Inv|"                          HEMCO_Diagn.rc

	# Turn @ into # characters for the benchmark simulation,
	# which should cause MAPL to skip reading these lines.
	# This is a workaround for a "input file to long" MAPL error.
	sed -i -e "s|@|#|"                               HISTORY.rc

	# Remove the first comment character on diagnostics
        sed -i -e "s|#'|'|"                              HISTORY.rc
    fi

    #------------------------------------------------------------------------
    # Standard settings
    #------------------------------------------------------------------------
    if [[ "x${sim_extra_option}" == "xnone" ]]; then
	sed -i -e 's/@//' HISTORY.rc
    fi

    #------------------------------------------------------------------------
    # Complex SOA settings
    #------------------------------------------------------------------------
    if [[ "x${sim_extra_option}" == "xbenchmark" ]] || \
       [[ ${sim_extra_option}    =~ "complexSOA" ]] || \
       [[ "x${sim_extra_option}" == "xAPM"       ]]; then

        # Turn on complex SOA option in input.geos
        replace_colon_sep_val "Online COMPLEX SOA" T input.geos

        # Add complex SOA species in input.geos
        prev_line="Species name            : ALK4"
        new_line="\Species name            : ASOA1\n\
Species name            : ASOA2\n\
Species name            : ASOA3\n\
Species name            : ASOAN\n\
Species name            : ASOG1\n\
Species name            : ASOG2\n\
Species name            : ASOG3"
        sed -i -e "/${prev_line}/a ${new_line}" input.geos

        prev_line="Species name            : TOLU"
        new_line="\Species name            : TSOA0\n\
Species name            : TSOA1\n\
Species name            : TSOA2\n\
Species name            : TSOA3\n\
Species name            : TSOG0\n\
Species name            : TSOG1\n\
Species name            : TSOG2\n\
Species name            : TSOG3"
        sed -i -e "/${prev_line}/a ${new_line}" input.geos

	sed -i -e 's/@//' HISTORY.rc
    fi

    # For complexSOA only, remove SOAP and SOAS species from input.geos
    if [[ ${sim_extra_option} =~ "complexSOA" ]]; then
	remove_text "Species name            : SOAP" input.geos
	remove_text "Species name            : SOAS" input.geos
    fi

    #------------------------------------------------------------------------
    # Semivolatile POA settings
    #------------------------------------------------------------------------
    if [[ "x${sim_extra_option}" == "xcomplexSOA_SVPOA" ]]; then

        # Turn on semivolatile POA option in input.geos
        replace_colon_sep_val "=> Semivolatile POA?" T input.geos

	# Remove non-SVPOA species from input.geos
        remove_text "Species name            : OCPI" input.geos
        remove_text "Species name            : OCPO" input.geos

        # Add semivolatile POA species in input.geos
        prev_line="Species name            : N2O5"
        new_line="\Species name            : NAP"
        sed -i -e "/${prev_line}/a ${new_line}" input.geos

        prev_line="Species name            : OIO"
        new_line="\Species name            : OPOA1\n\
Species name            : OPOA2\n\
Species name            : OPOG1\n\
Species name            : OPOG2"
        sed -i -e "/${prev_line}/a ${new_line}" input.geos

        prev_line="Species name            : PIP"
        new_line="\Species name            : POA1\n\
Species name            : POA2\n\
Species name            : POG1\n\
Species name            : POG2"
        sed -i -e "/${prev_line}/a ${new_line}" input.geos

	sed -i -e 's/@//' HISTORY.rc
    fi

    #------------------------------------------------------------------------
    # Acid uptake settings
    #------------------------------------------------------------------------
    if [[ "x${sim_extra_option}" == "xaciduptake" ]]; then
        replace_colon_sep_val "DustAlk" on HEMCO_Config.rc
        replace_colon_sep_val "=> Acidic uptake ?" T input.geos

        # Add DSTAL* species after DST4
        prev_line="Species name            : DST4"
        new_line="\Species name            : DSTAL1\n\
Species name            : DSTAL2\n\
Species name            : DSTAL3\n\
Species name            : DSTAL4"
        sed -i -e "/${prev_line}/a ${new_line}" input.geos

	# Add NITD* species after NITs.  NOTE: This is non-alphabetical,
	# but avoids double-adding these species after NIT and NITs.
        prev_line="Species name            : NITs"
        new_line="\Species name            : NITD1\n\
Species name            : NITD2\n\
Species name            : NITD3\n\
Species name            : NITD4"
        sed -i -e "/${prev_line}/a ${new_line}" input.geos

	# Add SO4* species after SO4s.  NOTE: This is non-alphabetical,
	# but avoids double-adding these species after SO4 and SO4s.
        prev_line="Species name            : SO4s"
        new_line="\Species name            : SO4D1\n\
Species name            : SO4D2\n\
Species name            : SO4D3\n\
Species name            : SO4D4"
        sed -i -e "/${prev_line}/a ${new_line}" input.geos

	sed -i -e 's/@//' HISTORY.rc
    fi

    #------------------------------------------------------------------------
    # Marine POA settings
    #------------------------------------------------------------------------
    if [[ "x${sim_extra_option}" == "xmarinePOA" ]]; then
        replace_colon_sep_val "SeaSalt"                 on HEMCO_Config.rc
        replace_colon_sep_val " => MARINE ORG AEROSOLS" T  input.geos

        # Add marine POA species to input.geos
        prev_line="Species name            : MONITU"
        new_line="\Species name            : MOPI\n\
Species name            : MOPO"
        sed -i -e "/${prev_line}/a ${new_line}" input.geos

	sed -i -e "s|@||"                      HISTORY.rc
    fi

    #------------------------------------------------------------------------
    # RRTMG settings
    #------------------------------------------------------------------------
    if [[ ${sim_extra_option} = "RRTMG" ]]; then

        replace_colon_sep_val "Turn on RRTMG?"       T input.geos
        replace_colon_sep_val "Calculate LW fluxes?" T input.geos
        replace_colon_sep_val "Calculate SW fluxes?" T input.geos
        replace_colon_sep_val "Clear-sky flux?"      T input.geos
        replace_colon_sep_val "All-sky flux?"        T input.geos
        replace_colon_sep_val "--> RRTMG"         true HEMCO_Config.rc

	sed -i -e "s|@||"                              HISTORY.rc
        sed -i -e "s|##'RRTMG'|'RRTMG'|"               HISTORY.rc
        printf "\nWARNING: All RRTMG run options are enabled which will significantly slow down the model!"
        printf "\nEdit input.geos and HISTORY.rc in your new run directory to customize options to only"
        printf "\nwhat you need.\n"
    fi

    #------------------------------------------------------------------------
    # TOMAS settings
    #------------------------------------------------------------------------
    if [[ ${sim_extra_option} =~ "TOMAS" ]]; then
        replace_colon_sep_val " => Use non-local PBL?"   F    input.geos
        replace_colon_sep_val "Use linear. strat. chem?" F    input.geos
        replace_colon_sep_val "=> Online O3 from model"  F    input.geos
        replace_colon_sep_val "TOMAS_Jeagle"             on   HEMCO_Config.rc
        # Add TOMAS species to input.geos
        prev_line="Species name            : XYLE"
        new_line="\Species name            : H2SO4\n\
Species name            : NK1\n\
Species name            : NK2\n\
Species name            : NK3\n\
Species name            : NK4\n\
Species name            : NK5\n\
Species name            : NK6\n\
Species name            : NK7\n\
Species name            : NK8\n\
Species name            : NK9\n\
Species name            : NK10\n\
Species name            : NK11\n\
Species name            : NK12\n\
Species name            : NK13\n\
Species name            : NK14\n\
Species name            : NK15\n\
Species name            : SF1\n\
Species name            : SF2\n\
Species name            : SF3\n\
Species name            : SF4\n\
Species name            : SF5\n\
Species name            : SF6\n\
Species name            : SF7\n\
Species name            : SF8\n\
Species name            : SF9\n\
Species name            : SF10\n\
Species name            : SF11\n\
Species name            : SF12\n\
Species name            : SF13\n\
Species name            : SF14\n\
Species name            : SF15\n\
Species name            : SS1\n\
Species name            : SS2\n\
Species name            : SS3\n\
Species name            : SS4\n\
Species name            : SS5\n\
Species name            : SS6\n\
Species name            : SS7\n\
Species name            : SS8\n\
Species name            : SS9\n\
Species name            : SS10\n\
Species name            : SS11\n\
Species name            : SS12\n\
Species name            : SS13\n\
Species name            : SS14\n\
Species name            : SS15\n\
Species name            : ECOB1\n\
Species name            : ECOB2\n\
Species name            : ECOB3\n\
Species name            : ECOB4\n\
Species name            : ECOB5\n\
Species name            : ECOB6\n\
Species name            : ECOB7\n\
Species name            : ECOB8\n\
Species name            : ECOB9\n\
Species name            : ECOB10\n\
Species name            : ECOB11\n\
Species name            : ECOB12\n\
Species name            : ECOB13\n\
Species name            : ECOB14\n\
Species name            : ECOB15\n\
Species name            : ECIL1\n\
Species name            : ECIL2\n\
Species name            : ECIL3\n\
Species name            : ECIL4\n\
Species name            : ECIL5\n\
Species name            : ECIL6\n\
Species name            : ECIL7\n\
Species name            : ECIL8\n\
Species name            : ECIL9\n\
Species name            : ECIL10\n\
Species name            : ECIL11\n\
Species name            : ECIL12\n\
Species name            : ECIL13\n\
Species name            : ECIL14\n\
Species name            : ECIL15\n\
Species name            : OCOB1\n\
Species name            : OCOB2\n\
Species name            : OCOB3\n\
Species name            : OCOB4\n\
Species name            : OCOB5\n\
Species name            : OCOB6\n\
Species name            : OCOB7\n\
Species name            : OCOB8\n\
Species name            : OCOB9\n\
Species name            : OCOB10\n\
Species name            : OCOB11\n\
Species name            : OCOB12\n\
Species name            : OCOB13\n\
Species name            : OCOB14\n\
Species name            : OCOB15\n\
Species name            : OCIL1\n\
Species name            : OCIL2\n\
Species name            : OCIL3\n\
Species name            : OCIL4\n\
Species name            : OCIL5\n\
Species name            : OCIL6\n\
Species name            : OCIL7\n\
Species name            : OCIL8\n\
Species name            : OCIL9\n\
Species name            : OCIL10\n\
Species name            : OCIL11\n\
Species name            : OCIL12\n\
Species name            : OCIL13\n\
Species name            : OCIL14\n\
Species name            : OCIL15\n\
Species name            : DUST1\n\
Species name            : DUST2\n\
Species name            : DUST3\n\
Species name            : DUST4\n\
Species name            : DUST5\n\
Species name            : DUST6\n\
Species name            : DUST7\n\
Species name            : DUST8\n\
Species name            : DUST9\n\
Species name            : DUST10\n\
Species name            : DUST11\n\
Species name            : DUST12\n\
Species name            : DUST13\n\
Species name            : DUST14\n\
Species name            : DUST15\n\
Species name            : AW1\n\
Species name            : AW2\n\
Species name            : AW3\n\
Species name            : AW4\n\
Species name            : AW5\n\
Species name            : AW6\n\
Species name            : AW7\n\
Species name            : AW8\n\
Species name            : AW9\n\
Species name            : AW10\n\
Species name            : AW11\n\
Species name            : AW12\n\
Species name            : AW13\n\
Species name            : AW14\n\
Species name            : AW15"
        sed -i -e "/${prev_line}/a ${new_line}" input.geos

        if [[ ${sim_extra_option} = "TOMAS40" ]]; then
    	prev_line="Species name            : NK15"
    	new_line="\Species name            : NK16\n\
Species name            : NK17\n\
Species name            : NK18\n\
Species name            : NK19\n\
Species name            : NK20\n\
Species name            : NK21\n\
Species name            : NK22\n\
Species name            : NK23\n\
Species name            : NK24\n\
Species name            : NK25\n\
Species name            : NK26\n\
Species name            : NK27\n\
Species name            : NK28\n\
Species name            : NK29\n\
Species name            : NK30\n\
Species name            : NK31\n\
Species name            : NK32\n\
Species name            : NK33\n\
Species name            : NK34\n\
Species name            : NK35\n\
Species name            : NK36\n\
Species name            : NK37\n\
Species name            : NK38\n\
Species name            : NK39\n\
Species name            : NK40"
    	sed -i -e "/${prev_line}/a ${new_line}" input.geos

	prev_line="Species name            : SF15"
    	new_line="\Species name            : SF16\n\
Species name            : SF17\n\
Species name            : SF18\n\
Species name            : SF19\n\
Species name            : SF20\n\
Species name            : SF21\n\
Species name            : SF22\n\
Species name            : SF23\n\
Species name            : SF24\n\
Species name            : SF25\n\
Species name            : SF26\n\
Species name            : SF27\n\
Species name            : SF28\n\
Species name            : SF29\n\
Species name            : SF30\n\
Species name            : SF31\n\
Species name            : SF32\n\
Species name            : SF33\n\
Species name            : SF34\n\
Species name            : SF35\n\
Species name            : SF36\n\
Species name            : SF37\n\
Species name            : SF38\n\
Species name            : SF39\n\
Species name            : SF40"
    	sed -i -e "/${prev_line}/a ${new_line}" input.geos

    	prev_line="Species name            : SS15"
    	new_line="\Species name            : SS16\n\
Species name            : SS17\n\
Species name            : SS18\n\
Species name            : SS19\n\
Species name            : SS20\n\
Species name            : SS21\n\
Species name            : SS22\n\
Species name            : SS23\n\
Species name            : SS24\n\
Species name            : SS25\n\
Species name            : SS26\n\
Species name            : SS27\n\
Species name            : SS28\n\
Species name            : SS29\n\
Species name            : SS30\n\
Species name            : SS31\n\
Species name            : SS32\n\
Species name            : SS33\n\
Species name            : SS34\n\
Species name            : SS35\n\
Species name            : SS36\n\
Species name            : SS37\n\
Species name            : SS38\n\
Species name            : SS39\n\
Species name            : SS40"
    	sed -i -e "/${prev_line}/a ${new_line}" input.geos

    	prev_line="Species name            : ECOB15"
    	new_line="\Species name            : ECOB16\n\
Species name            : ECOB17\n\
Species name            : ECOB18\n\
Species name            : ECOB19\n\
Species name            : ECOB20\n\
Species name            : ECOB21\n\
Species name            : ECOB22\n\
Species name            : ECOB23\n\
Species name            : ECOB24\n\
Species name            : ECOB25\n\
Species name            : ECOB26\n\
Species name            : ECOB27\n\
Species name            : ECOB28\n\
Species name            : ECOB29\n\
Species name            : ECOB30\n\
Species name            : ECOB31\n\
Species name            : ECOB32\n\
Species name            : ECOB33\n\
Species name            : ECOB34\n\
Species name            : ECOB35\n\
Species name            : ECOB36\n\
Species name            : ECOB37\n\
Species name            : ECOB38\n\
Species name            : ECOB39\n\
Species name            : ECOB40"
    	sed -i -e "/${prev_line}/a ${new_line}" input.geos

    	prev_line="Species name            : ECIL15"
    	new_line="\Species name            : ECIL16\n\
Species name            : ECIL17\n\
Species name            : ECIL18\n\
Species name            : ECIL19\n\
Species name            : ECIL20\n\
Species name            : ECIL21\n\
Species name            : ECIL22\n\
Species name            : ECIL23\n\
Species name            : ECIL24\n\
Species name            : ECIL25\n\
Species name            : ECIL26\n\
Species name            : ECIL27\n\
Species name            : ECIL28\n\
Species name            : ECIL29\n\
Species name            : ECIL30\n\
Species name            : ECIL31\n\
Species name            : ECIL32\n\
Species name            : ECIL33\n\
Species name            : ECIL34\n\
Species name            : ECIL35\n\
Species name            : ECIL36\n\
Species name            : ECIL37\n\
Species name            : ECIL38\n\
Species name            : ECIL39\n\
Species name            : ECIL40"
    	sed -i -e "/${prev_line}/a ${new_line}" input.geos

    	prev_line="Species name            : OCOB15"
    	new_line="\Species name            : OCOB16\n\
Species name            : OCOB17\n\
Species name            : OCOB18\n\
Species name            : OCOB19\n\
Species name            : OCOB20\n\
Species name            : OCOB21\n\
Species name            : OCOB22\n\
Species name            : OCOB23\n\
Species name            : OCOB24\n\
Species name            : OCOB25\n\
Species name            : OCOB26\n\
Species name            : OCOB27\n\
Species name            : OCOB28\n\
Species name            : OCOB29\n\
Species name            : OCOB30\n\
Species name            : OCOB31\n\
Species name            : OCOB32\n\
Species name            : OCOB33\n\
Species name            : OCOB34\n\
Species name            : OCOB35\n\
Species name            : OCOB36\n\
Species name            : OCOB37\n\
Species name            : OCOB38\n\
Species name            : OCOB39\n\
Species name            : OCOB40"
    	sed -i -e "/${prev_line}/a ${new_line}" input.geos

    	prev_line="Species name            : OCIL15"
    	new_line="\Species name            : OCIL16\n\
Species name            : OCIL17\n\
Species name            : OCIL18\n\
Species name            : OCIL19\n\
Species name            : OCIL20\n\
Species name            : OCIL21\n\
Species name            : OCIL22\n\
Species name            : OCIL23\n\
Species name            : OCIL24\n\
Species name            : OCIL25\n\
Species name            : OCIL26\n\
Species name            : OCIL27\n\
Species name            : OCIL28\n\
Species name            : OCIL29\n\
Species name            : OCIL30\n\
Species name            : OCIL31\n\
Species name            : OCIL32\n\
Species name            : OCIL33\n\
Species name            : OCIL34\n\
Species name            : OCIL35\n\
Species name            : OCIL36\n\
Species name            : OCIL37\n\
Species name            : OCIL38\n\
Species name            : OCIL39\n\
Species name            : OCIL40"
    	sed -i -e "/${prev_line}/a ${new_line}" input.geos

    	prev_line="Species name            : DUST15"
    	new_line="\Species name            : DUST16\n\
Species name            : DUST17\n\
Species name            : DUST18\n\
Species name            : DUST19\n\
Species name            : DUST20\n\
Species name            : DUST21\n\
Species name            : DUST22\n\
Species name            : DUST23\n\
Species name            : DUST24\n\
Species name            : DUST25\n\
Species name            : DUST26\n\
Species name            : DUST27\n\
Species name            : DUST28\n\
Species name            : DUST29\n\
Species name            : DUST30\n\
Species name            : DUST31\n\
Species name            : DUST32\n\
Species name            : DUST33\n\
Species name            : DUST34\n\
Species name            : DUST35\n\
Species name            : DUST36\n\
Species name            : DUST37\n\
Species name            : DUST38\n\
Species name            : DUST39\n\
Species name            : DUST40"
    	sed -i -e "/${prev_line}/a ${new_line}" input.geos

    	prev_line="Species name            : AW15"
    	new_line="\Species name            : AW16\n\
Species name            : AW17\n\
Species name            : AW18\n\
Species name            : AW19\n\
Species name            : AW20\n\
Species name            : AW21\n\
Species name            : AW22\n\
Species name            : AW23\n\
Species name            : AW24\n\
Species name            : AW25\n\
Species name            : AW26\n\
Species name            : AW27\n\
Species name            : AW28\n\
Species name            : AW29\n\
Species name            : AW30\n\
Species name            : AW31\n\
Species name            : AW32\n\
Species name            : AW33\n\
Species name            : AW34\n\
Species name            : AW35\n\
Species name            : AW36\n\
Species name            : AW37\n\
Species name            : AW38\n\
Species name            : AW39\n\
Species name            : AW40"
    	sed -i -e "/${prev_line}/a ${new_line}" input.geos
	sed -i -e 's/@//' HISTORY.rc
        fi
    fi

    if [[ ${sim_extra_option} = "APM" ]]; then

        # Add APM species to input.geos
        prev_line="Species name            : XYLE"
        new_line="\Species name            : APMBCBIN01\n\
Species name            : APMBCBIN02\n\
Species name            : APMBCBIN03\n\
Species name            : APMBCBIN04\n\
Species name            : APMBCBIN05\n\
Species name            : APMBCBIN06\n\
Species name            : APMBCBIN07\n\
Species name            : APMBCBIN08\n\
Species name            : APMBCBIN09\n\
Species name            : APMBCBIN10\n\
Species name            : APMBCBIN11\n\
Species name            : APMBCBIN12\n\
Species name            : APMBCBIN13\n\
Species name            : APMBCBIN14\n\
Species name            : APMBCBIN15\n\
Species name            : APMCTBC1\n\
Species name            : APMCTBC2\n\
Species name            : APMCTDST1\n\
Species name            : APMCTDST2\n\
Species name            : APMCTOC1\n\
Species name            : APMCTOC2\n\
Species name            : APMCTSEA1\n\
Species name            : APMCTSEA2\n\
Species name            : APMDSTBIN01\n\
Species name            : APMDSTBIN02\n\
Species name            : APMDSTBIN03\n\
Species name            : APMDSTBIN04\n\
Species name            : APMDSTBIN05\n\
Species name            : APMDSTBIN06\n\
Species name            : APMDSTBIN07\n\
Species name            : APMDSTBIN08\n\
Species name            : APMDSTBIN09\n\
Species name            : APMDSTBIN10\n\
Species name            : APMDSTBIN11\n\
Species name            : APMDSTBIN12\n\
Species name            : APMDSTBIN13\n\
Species name            : APMDSTBIN14\n\
Species name            : APMDSTBIN15\n\
Species name            : APMH2SO4\n\
Species name            : APMLVSOA\n\
Species name            : APMLVSOG\n\
Species name            : APMOCBIN01\n\
Species name            : APMOCBIN02\n\
Species name            : APMOCBIN03\n\
Species name            : APMOCBIN04\n\
Species name            : APMOCBIN05\n\
Species name            : APMOCBIN06\n\
Species name            : APMOCBIN07\n\
Species name            : APMOCBIN08\n\
Species name            : APMOCBIN09\n\
Species name            : APMOCBIN10\n\
Species name            : APMOCBIN11\n\
Species name            : APMOCBIN12\n\
Species name            : APMOCBIN13\n\
Species name            : APMOCBIN14\n\
Species name            : APMOCBIN15\n\
Species name            : APMSEABIN01\n\
Species name            : APMSEABIN02\n\
Species name            : APMSEABIN03\n\
Species name            : APMSEABIN04\n\
Species name            : APMSEABIN05\n\
Species name            : APMSEABIN06\n\
Species name            : APMSEABIN07\n\
Species name            : APMSEABIN08\n\
Species name            : APMSEABIN09\n\
Species name            : APMSEABIN10\n\
Species name            : APMSEABIN11\n\
Species name            : APMSEABIN12\n\
Species name            : APMSEABIN13\n\
Species name            : APMSEABIN14\n\
Species name            : APMSEABIN15\n\
Species name            : APMSEABIN16\n\
Species name            : APMSEABIN17\n\
Species name            : APMSEABIN18\n\
Species name            : APMSEABIN19\n\
Species name            : APMSEABIN20\n\
Species name            : APMSPBIN01\n\
Species name            : APMSPBIN02\n\
Species name            : APMSPBIN03\n\
Species name            : APMSPBIN04\n\
Species name            : APMSPBIN05\n\
Species name            : APMSPBIN06\n\
Species name            : APMSPBIN07\n\
Species name            : APMSPBIN08\n\
Species name            : APMSPBIN09\n\
Species name            : APMSPBIN10\n\
Species name            : APMSPBIN11\n\
Species name            : APMSPBIN12\n\
Species name            : APMSPBIN13\n\
Species name            : APMSPBIN14\n\
Species name            : APMSPBIN15\n\
Species name            : APMSPBIN16\n\
Species name            : APMSPBIN17\n\
Species name            : APMSPBIN18\n\
Species name            : APMSPBIN19\n\
Species name            : APMSPBIN20\n\
Species name            : APMSPBIN21\n\
Species name            : APMSPBIN22\n\
Species name            : APMSPBIN23\n\
Species name            : APMSPBIN24\n\
Species name            : APMSPBIN25\n\
Species name            : APMSPBIN26\n\
Species name            : APMSPBIN27\n\
Species name            : APMSPBIN28\n\
Species name            : APMSPBIN29\n\
Species name            : APMSPBIN30\n\
Species name            : APMSPBIN31\n\
Species name            : APMSPBIN32\n\
Species name            : APMSPBIN33\n\
Species name            : APMSPBIN34\n\
Species name            : APMSPBIN35\n\
Species name            : APMSPBIN36\n\
Species name            : APMSPBIN37\n\
Species name            : APMSPBIN38\n\
Species name            : APMSPBIN39\n\
Species name            : APMSPBIN40"
        sed -i -e "/${prev_line}/a ${new_line}" input.geos
	sed -i -e 's/@//' HISTORY.rc
    fi
}
#EOC
