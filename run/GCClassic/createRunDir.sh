#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: createRunDir.sh
#
# !DESCRIPTION: Creates a GEOS-Chem Classic run directory.
#\\
#\\
# !CALLING SEQUENCE:
#  ./createRunDir.sh [rundirname]
#
# !REMARKS:
#  If optional run directory name argument is not passed then the user
#  will be prompted to enter a name interactively, or choose to use the
#  default name gc_{grid_display}_{met}_{sim_name}_{sim_extra_option}.
#
# !REVISION HISTORY:
#  Initial version: M. Sulprizio, 6/24/2020 (based off GCHP/createRunDir.sh)
#  See the subsequent Git history with the gitk browser!
#------------------------------------------------------------------------------
#BOC

# Directory with GEOS-Chem Classic rundir scripts (i.e. this dir)
# Current directory
srcrundir=$(pwd -P)
cd ${srcrundir}

# GEOS-Chem "science codebase" directory
cd ../..
gcdir=$(pwd -P)

# GCClassic "wrapper" directory
cd ../../
wrapperdir=$(pwd -P)

# Return to directory w/ GEOS-Chem Classic rundir scripts
cd ${srcrundir}

# Source common bash functions from scripts in the run/shared folder
. ${gcdir}/run/shared/setupConfigFiles.sh      # Config file editing
. ${gcdir}/run/GCClassic/setupForRestarts.sh   # Functions for restart files
. ${gcdir}/run/shared/newUserRegistration.sh   # 1st-time user registration
. ${gcdir}/run/shared/singleCarbonSpecies.sh   # Single carbon species setup

# Initialize run directory variables
RUNDIR_VARS=""
RUNDIR_VARS+="RUNDIR_GC_MODE='GCClassic'\n"

# Define separator lines
thickline="\n===========================================================\n"
thinline="\n-----------------------------------------------------------\n"

printf "${thickline}GEOS-CHEM RUN DIRECTORY CREATION${thickline}"

#-----------------------------------------------------------------
# Export data root path in ~/.geoschem/config if file exists
#-----------------------------------------------------------------
if [[ -f ${HOME}/.geoschem/config ]]; then
    source ${HOME}/.geoschem/config
    if [[ ! -d ${GC_DATA_ROOT} ]]; then
	printf "\nWarning: Default root data directory does not exist!"
        printf "\nSet new path below or manually edit ${HOME}/.geoschem/config.\n"
    fi
else
    printf "${thinline}Define path to ExtData."
    printf "\nThis will be stored in ${HOME}/.geoschem/config for future automatic use.${thinline}"
    mkdir -p ${HOME}/.geoschem
fi

#-----------------------------------------------------------------
# One-time configuration of data root path in ~/.geoschem/config
#-----------------------------------------------------------------
if [[ -z "${GC_DATA_ROOT}" ]]; then
    printf "${thinline}Enter path for ExtData:${thinline}"
    valid_path=0
    while [ "$valid_path" -eq 0 ]; do
	read -e -p "${USER_PROMPT}" extdata
	if [[ ${extdata} = "q" ]]; then
	    printf "\nExiting.\n"
	    exit 1
	elif [[ ! -d ${extdata} ]]; then
            printf "\nERROR: ${extdata} does not exist. Enter a new path or hit q to quit.\n"
	else
	    valid_path=1
	    echo "export GC_DATA_ROOT=${extdata}" >> ${HOME}/.geoschem/config
            source ${HOME}/.geoschem/config
	fi
    done
fi

RUNDIR_VARS+="RUNDIR_DATA_ROOT=$GC_DATA_ROOT\n"

# --------------------------------------------------------------
# registration for first time users
# --------------------------------------------------------------
[[ -z "${GC_USER_REGISTERED}" ]] && registerNewUser "gcc"

#-----------------------------------------------------------------
# Ask user to select simulation type
#-----------------------------------------------------------------
printf "${thinline}Choose simulation type:${thinline}"
printf "   1. Full chemistry\n"
printf "   2. Aerosols only\n"
printf "   3. Carbon\n"
printf "   4. Hg\n"
printf "   5. POPs\n"
printf "   6. Tagged O3\n"
printf "   7. TransportTracers\n"
printf "   8. Trace metals\n"
printf "   9. CH4\n"
printf "  10. CO2\n"
printf "  11. Tagged CO\n"
valid_sim=0
while [ "${valid_sim}" -eq 0 ]; do
    read -p "${USER_PROMPT}" sim_num
    valid_sim=1
    if [[ ${sim_num} = "1" ]]; then
	sim_name=fullchem
    elif [[ ${sim_num} = "2" ]]; then
	sim_name=aerosol
    elif [[ ${sim_num} = "3" ]]; then
	sim_name=carbon
    elif [[ ${sim_num} = "4" ]]; then
	sim_name=Hg
    elif [[ ${sim_num} = "5" ]]; then
	sim_name=POPs
    elif [[ ${sim_num} = "6" ]]; then
	sim_name=tagO3
    elif [[ ${sim_num} = "7" ]]; then
	sim_name=TransportTracers
    elif [[ ${sim_num} = "8" ]]; then
	sim_name=metals
    elif [[ ${sim_num} = "9" ]]; then
	sim_name=CH4
    elif [[ ${sim_num} = "10" ]]; then
	sim_name=CO2
    elif [[ ${sim_num} = "11" ]]; then
	sim_name=tagCO
    else
        valid_sim=0
	printf "Invalid simulation option. Try again.\n"
    fi
done

RUNDIR_VARS+="RUNDIR_SIM_NAME=$sim_name\n"

#-----------------------------------------------------------------
# Ask user to specify full-chemistry simulation options
#-----------------------------------------------------------------
sim_extra_option=none

# Ask user to specify full chemistry simulation options
if [[ ${sim_name} = "fullchem" ]]; then

    printf "${thinline}Choose additional simulation option:${thinline}"
    printf "  1. Standard\n"
    printf "  2. Benchmark\n"
    printf "  3. Complex SOA\n"
    printf "  4. Marine POA\n"
    printf "  5. Acid uptake on dust\n"
    printf "  6. TOMAS\n"
    printf "  7. APM\n"
    printf "  8. RRTMG\n"
    valid_sim_option=0
    while [ "${valid_sim_option}" -eq 0 ]; do
	read -p "${USER_PROMPT}" sim_option
	valid_sim_option=1
	if [[ ${sim_option} = "1" ]]; then
	    sim_extra_option=none
	elif [[ ${sim_option} = "2" ]]; then
	    sim_extra_option="benchmark"
	elif [[ ${sim_option} = "3" ]]; then
	    printf "${thinline}Choose complex SOA option:${thinline}"
	    printf "  1. Complex SOA\n"
	    printf "  2. Complex SOA with semivolatile POA\n"
	    valid_soa=0
	    while [ "${valid_soa}" -eq 0 ]; do
		read -p "${USER_PROMPT}" soa_option
		valid_soa=1
		if [[ ${soa_option} = "1" ]]; then
		    sim_extra_option="complexSOA"
		elif [[ ${soa_option} = "2" ]]; then
		    sim_extra_option="complexSOA_SVPOA"
		else
		    valid_soa=0
		    printf "Invalid complex SOA option.Try again.\n"
		fi
	    done
	elif [[ ${sim_option} = "4" ]]; then
	    sim_extra_option="marinePOA"
	elif [[ ${sim_option} = "5" ]]; then
	    sim_extra_option="aciduptake"
	elif [[ ${sim_option} = "6" ]]; then
	    printf "${thinline}Choose TOMAS option:${thinline}"
	    printf "  1. TOMAS with 15 bins\n"
	    printf "  2. TOMAS with 40 bins\n"
	    valid_tomas=0
	    while [ "${valid_tomas}" -eq 0 ]; do
		read -p "${USER_PROMPT}" tomas_option
		valid_tomas=1
		if [[ ${tomas_option} = "1" ]]; then
		    sim_extra_option="TOMAS15"
		elif [[ ${tomas_option} = "2" ]]; then
		    sim_extra_option="TOMAS40"
		else
		    valid_tomas=0
		    printf "Invalid TOMAS option. Try again.\n"
		fi
	    done
	elif [[ ${sim_option} = "7" ]]; then
	    sim_extra_option="APM"
	elif [[ ${sim_option} = "8" ]]; then
	    sim_extra_option="RRTMG"
	else
	    valid_sim_option=0
	    printf "Invalid simulation option. Try again.\n"
	fi
    done

# Currently no transport tracer extra options
elif [[ ${sim_name} = "TransportTracers" ]]; then
    sim_extra_option=none

# Ask user to specify POPs simulation options
elif [[ ${sim_name} = "POPs" ]]; then
    printf "${thinline}Choose POPs type:${thinline}"
    printf "  1. BaP\n"
    printf "  2. PHE\n"
    printf "  3. PYR\n"
    valid_pops=0
    while [ "${valid_pops}" -eq 0 ]; do
	read -p "${USER_PROMPT}" pops_num
	valid_pops=1
	if [[ ${pops_num} = "1" ]]; then
	    sim_extra_option="BaP"
	elif [[ ${pops_num} = "2" ]]; then
	    sim_extra_option="PHE"
	elif [[ ${pops_num} = "3" ]]; then
	    sim_extra_option="PYR"
	else
	    valid_pops=0
	    printf "Invalid POPs type. Try again.\n"
	fi
    done

# Ask user to specify carbon simulation options
elif [[ "x${sim_name}" == "xcarbon" ]]; then
    printf "${thinline}Do you wish to use a single advected species?${thinline}"
    printf "  1. Use all species\n"
    printf "  2. Use CH4 only\n"
    printf "  3. Use CO2 only\n"
    printf "  4. Use CO only\n"
    printf "  5. Use OCS only\n"
    valid=0
    while [ "${valid}" -eq 0 ]; do
	read -p "${USER_PROMPT}" prompt
	valid=1
	if [[ "x${prompt}" == "x1" ]]; then
	    sim_extra_option="none"
	elif [[ "x${prompt}" == "x2" ]]; then
	    sim_extra_option="CH4"
	elif [[ "x${prompt}" == "x3" ]]; then
	    sim_extra_option="CO2"
	elif [[ "x${prompt}" == "x4" ]]; then
	    sim_extra_option="CO"
	elif [[ "x${prompt}" == "x5" ]]; then
	    sim_extra_option="OCS"
	else
	    valid=0
	    printf "Invalid selection. Try again.\n"
	fi
    done
fi

RUNDIR_VARS+="RUNDIR_SIM_EXTRA_OPTION=$sim_extra_option\n"

# Determine settings based on simulation type
SettingsDir="${gcdir}/run/shared/settings"
if [[ ${sim_extra_option} == "BaP" ]]; then
    RUNDIR_VARS+="$(cat ${SettingsDir}/POPs_BaP.txt)\n"
elif [[ ${sim_extra_option} == "PHE" ]]; then
    RUNDIR_VARS+="$(cat ${SettingsDir}/POPs_PHE.txt)\n"
elif [[ ${sim_extra_option} == "PYR" ]]; then
    RUNDIR_VARS+="$(cat ${SettingsDir}/POPs_PYR.txt)\n"
fi

if [[ ${sim_extra_option} == "benchmark"  ]] || \
   [[ ${sim_extra_option} =~ "complexSOA" ]] || \
   [[ ${sim_extra_option} == "APM"        ]]; then
    RUNDIR_VARS+="RUNDIR_COMPLEX_SOA='true '\n"
    if [[ ${sim_extra_option} == "complexSOA_SVPOA" ]]; then
	RUNDIR_VARS+="RUNDIR_SVPOA='true '\n"
    else
	RUNDIR_VARS+="RUNDIR_SVPOA='false'\n"
    fi
else
    RUNDIR_VARS+="RUNDIR_COMPLEX_SOA='false'\n"
    RUNDIR_VARS+="RUNDIR_SVPOA='false'\n"
fi

if [[ ${sim_extra_option} == "aciduptake" ]]; then
    RUNDIR_VARS+="RUNDIR_DUSTALK_EXT='on '\n"
    RUNDIR_VARS+="RUNDIR_ACID_UPTAKE='true '\n"
else
    RUNDIR_VARS+="RUNDIR_DUSTALK_EXT='off'\n"
    RUNDIR_VARS+="RUNDIR_ACID_UPTAKE='false'\n"
fi

if [[ ${sim_extra_option} == "marinePOA" ]]; then
    RUNDIR_VARS+="RUNDIR_MARINE_POA='true '\n"
else
    RUNDIR_VARS+="RUNDIR_MARINE_POA='false'\n"
fi

if [[ ${sim_extra_option} == "RRTMG" ]]; then
    RUNDIR_VARS+="RUNDIR_RRTMG_OPTS='true '\n"
    RUNDIR_VARS+="RUNDIR_USE_RRTMG='true '\n"
else
    RUNDIR_VARS+="RUNDIR_RRTMG_OPTS='false'\n"
    RUNDIR_VARS+="RUNDIR_USE_RRTMG='false'\n"
fi

if [[ ${sim_extra_option} =~ "TOMAS" ]]; then
    RUNDIR_VARS+="RUNDIR_USE_NLPBL='false'\n"
    RUNDIR_VARS+="RUNDIR_USE_ONLINE_O3='false'\n"
else
    RUNDIR_VARS+="RUNDIR_USE_NLPBL='true '\n"
    RUNDIR_VARS+="RUNDIR_USE_ONLINE_O3='true '\n"
fi

#-----------------------------------------------------------------
# Ask user to select meteorology source
#-----------------------------------------------------------------
printf "${thinline}Choose meteorology source:${thinline}"
printf "  1. MERRA-2 (Recommended)\n"
printf "  2. GEOS-FP \n"
printf "  3. GEOS-IT (Beta release)\n"
printf "  4. GISS ModelE2.1 (GCAP 2.0)\n"

valid_met=0
while [ "${valid_met}" -eq 0 ]; do
    read -p "${USER_PROMPT}" met_num
    valid_met=1
    if [[ ${met_num} = "1" ]]; then
	met="merra2"
	shared_met_settings=${gcdir}/run/shared/settings/merra2.txt
	RUNDIR_VARS+="RUNDIR_MET_FIELD_CONFIG='HEMCO_Config.rc.gmao_metfields'\n"
    elif [[ ${met_num} = "2" ]]; then
	met="geosfp"
	shared_met_settings=${gcdir}/run/shared/settings/geosfp/geosfp.preprocessed_ll.txt
	RUNDIR_VARS+="RUNDIR_MET_FIELD_CONFIG='HEMCO_Config.rc.gmao_metfields'\n"

	# Print warning about GEOS-FP and require user to acknowledge it.
	fp_msg="WARNING: The convection scheme used to generate archived GEOS-FP meteorology \nfiles changed from RAS to Grell-Freitas starting June 1 2020 with impact on \nvertical transport. Discussion and analysis of the impact is available at \ngithub.com/geoschem/geos-chem/issues/1409. In addition, there is a bug in \nconvective precipitation flux following the switch where all values are zero \nin the input files. This bug is addressed by computing fluxes online for runs \nstarting on or after June 1 2020. The fix does not extend to the case of running \nacross the time boundary. Due to these issues we recommend splitting up GEOS-FP \nruns in time such that a single simulation does not span the switch. Configure \none run to end on June 1 2020 and then use its output restart to start another \nrun on June 1. Alternatively consider using MERRA2. If you wish to use a \nGEOS-FP meteorology year different from your simulation year please create a \nGEOS-Chem GitHub issue for assistance to avoid accidentally using zero \nconvective precipitation flux.\n"
	printf "\n${fp_msg}\n"
	printf "This warning will be printed to run directory file warnings.txt.\n"
	printf "${thinline}Enter y to acknowledge and proceed, or q to quit:${thinline}"
	valid_fp_accept=0
	while [ "${valid_fp_accept}" -eq 0 ]; do
	    read -p "${USER_PROMPT}" fp_accept
	    valid_fp_accept=1
	    if [[ ${fp_accept} = "y" ]]; then
		x=0
	    elif [[ ${fp_accept} = "q" ]]; then
		exit
	    else
		valid_fp_accept=0
		printf "Invalid option. Try again.\n"
	    fi
	done

    elif [[ ${met_num} = "3" ]]; then
	met="geosit"
	shared_met_settings=${gcdir}/run/shared/settings/geosit/geosit.preprocessed_ll.txt
	RUNDIR_VARS+="RUNDIR_MET_FIELD_CONFIG='HEMCO_Config.rc.gmao_metfields'\n"
    elif [[ ${met_num} = "4" ]]; then
	met="ModelE2.1"
	shared_met_settings=${gcdir}/run/shared/settings/modele2.1.txt
	RUNDIR_VARS+="RUNDIR_MET_FIELD_CONFIG='HEMCO_Config.rc.gcap2_metfields'\n"
    else
	valid_met=0
	printf "Invalid meteorology option. Try again.\n"
    fi
done

if [[ ${met} = "ModelE2.1" ]]; then
    printf "${thinline}Choose scenario (presently available years in parentheses):${thinline}"
    printf "  1. Historical (1851-1860; 2001-2014)\n"
    printf "  2. Historical nudged to MERRA-2 (2001-2014)\n"
    printf "  3. SSP1-1.9 (2040-2049; 2090-2099)\n"
    printf "  4. SSP1-2.6 (2040-2049; 2090-2099)\n"
    printf "  5. SSP4-3.4 (2040-2049; 2090-2099)\n"
    printf "  6. SSP2-4.5 (2040-2049; 2090-2099)\n"
    printf "  7. SSP4-6.0 (2040-2049; 2090-2099)\n"
    printf "  8. SSP3-7.0 (2040-2049; 2090-2099)\n"
    printf "  9. SSP5-8.5 (2040-2049; 2090-2099)\n"

    valid_scen=0
    while [ "${valid_scen}" -eq 0 ]; do
	read -p "${USER_PROMPT}" scen_num
	valid_scen=1
	if [[ ${scen_num} = "1" ]]; then
	    scenario="HIST"
	    runid="E213f10aF40oQ40"
            gissres="F40"
	    vertres="40L"
	    RUNDIR_VARS+="RUNDIR_MET_AVAIL='# 1851-1860; 2001-2014'\n"
	elif [[ ${scen_num} = "2" ]]; then
	    scenario="HIST"
	    runid="E213f10aF40oQ40nudge"
            gissres="F40"
	    vertres="40L"
	    RUNDIR_VARS+="RUNDIR_MET_AVAIL='# 2001-2014'\n"
	elif [[ ${scen_num} = "3" ]]; then
	    scenario="SSP119"
	    runid="E213SSP119aF40oQ40"
            gissres="F40"
	    vertres="40L"
	    RUNDIR_VARS+="RUNDIR_MET_AVAIL='# 2040-2049; 2090-2099'\n"
	elif [[ ${scen_num} = "4" ]]; then
	    scenario="SSP126"
	    runid="E213SSP126aF40oQ40"
            gissres="F40"
	    vertres="40L"
	    RUNDIR_VARS+="RUNDIR_MET_AVAIL='# 2040-2049; 2090-2099'\n"
	elif [[ ${scen_num} = "5" ]]; then
	    scenario="SSP434"
	    runid="E213SSP434aF40oQ40"
            gissres="F40"
	    vertres="40L"
	    RUNDIR_VARS+="RUNDIR_MET_AVAIL='# 2040-2049; 2090-2099'\n"
	elif [[ ${scen_num} = "6" ]]; then
	    scenario="SSP245"
	    runid="E213SSP245aF40oQ40"
            gissres="F40"
	    vertres="40L"
	    RUNDIR_VARS+="RUNDIR_MET_AVAIL='# 2040-2049; 2090-2099'\n"
	elif [[ ${scen_num} = "7" ]]; then
	    scenario="SSP460"
	    runid="E213SSP460aF40oQ40"
            gissres="F40"
	    vertres="40L"
	    RUNDIR_VARS+="RUNDIR_MET_AVAIL='# 2040-2049; 2090-2099'\n"
	elif [[ ${scen_num} = "8" ]]; then
	    scenario="SSP370"
	    runid="E213SSP370aF40oQ40"
            gissres="F40"
	    vertres="40L"
	    RUNDIR_VARS+="RUNDIR_MET_AVAIL='# 2040-2049; 2090-2099'\n"
	elif [[ ${scen_num} = "9" ]]; then
	    scenario="SSP585"
	    runid="E213SSP585aF40oQ40"
            gissres="F40"
	    vertres="40L"
	    RUNDIR_VARS+="RUNDIR_MET_AVAIL='# 2040-2049; 2090-2099'\n"
	else
  	    valid_scen=0
	    printf "Invalid GCAP 2.0 scenario. Try again.\n"
	fi
	RUNDIR_VARS+="RUNDIR_GCAP2_SCENARIO='$scenario'\n"
        RUNDIR_VARS+="RUNDIR_GISS_RES='$gissres'\n"
	RUNDIR_VARS+="RUNDIR_GCAP2_VERTRES='$vertres'\n"
	RUNDIR_VARS+="RUNDIR_GCAP2_RUNID='$runid'\n"
    done

    RUNDIR_VARS+="RUNDIR_USE_AEIC='false'\n"

    # Define the volcano paths for the HEMCO_Config.rc file
    # NOTE: Benchmark simulations always use the climatological emissions!
    if [[ "x${sim_name}" == "xfullchem" ]]  ||  \
       [[ "x${sim_name}" == "xaerosol"  ]]; then
        RUNDIR_VARS+="RUNDIR_VOLC_CLIMATOLOGY='\$ROOT/VOLCANO/v2024-04/so2_volcanic_emissions_CARN_v202401.degassing_only.rc'\n"

	if [[ "x${sim_extra_option}" == "xbenchmark" ]]; then
	    RUNDIR_VARS+="RUNDIR_VOLC_TABLE='\$ROOT/VOLCANO/v2024-04/so2_volcanic_emissions_CARN_v202401.degassing_only.rc'\n"
	else
	    RUNDIR_VARS+="RUNDIR_VOLC_TABLE='\$ROOT/VOLCANO/v2024-04/\$YYYY/\$MM/so2_volcanic_emissions_Carns.\$YYYY\$MM\$DD.rc'\n"
	fi
    fi

else
    RUNDIR_VARS+="RUNDIR_GCAP2_SCENARIO='not_used'\n"
    RUNDIR_VARS+="RUNDIR_GISS_RES='not_used'\n"
    RUNDIR_VARS+="RUNDIR_GCAP2_VERTRES='not_used'\n"
    RUNDIR_VARS+="RUNDIR_GCAP2_RUNID='not_used'\n"
    RUNDIR_VARS+="RUNDIR_MET_AVAIL='# 1980-2021'\n"

    RUNDIR_VARS+="RUNDIR_USE_AEIC='true'\n"

    # Define the volcano paths for the HEMCO_Config.rc file
    # NOTE: Benchmark simulations always use the climatological emissions!
    if [[ "x${sim_name}" == "xfullchem" ]]  ||  \
       [[ "x${sim_name}" == "xaerosol"  ]]; then
	RUNDIR_VARS+="RUNDIR_VOLC_CLIMATOLOGY='\$ROOT/VOLCANO/v2024-04/so2_volcanic_emissions_CARN_v202401.degassing_only.rc'\n"

	if [[ "x${sim_extra_option}" == "xbenchmark" ]]; then
	    RUNDIR_VARS+="RUNDIR_VOLC_TABLE='\$ROOT/VOLCANO/v2024-04/so2_volcanic_emissions_CARN_v202401.degassing_only.rc'\n"
	else
	    RUNDIR_VARS+="RUNDIR_VOLC_TABLE='\$ROOT/VOLCANO/v2024-04/\$YYYY/\$MM/so2_volcanic_emissions_Carns.\$YYYY\$MM\$DD.rc'\n"
	fi
    fi

 fi

# Turn off MEGAN for the aerosol-only simulation
if [[ "x${sim_name}" == "xaerosol"  ]]; then
    RUNDIR_VARS+="RUNDIR_MEGAN_EXT='off'\n"
fi

#-----------------------------------------------------------------
# Ask user to select horizontal resolution
#-----------------------------------------------------------------
printf "${thinline}Choose horizontal resolution:${thinline}"
if [[ "x${met}" == "xModelE2.1" || "x${met}" == "xModelE2.2" ]]; then
    printf "  1. 4.0  x 5.0 *\n"
    printf "  2. 2.0  x 2.5\n"
    printf "  3. 0.5  x 0.625 *\n"
    printf "  4. 0.25 x 0.3125 *${thinline}"
    printf "  \n* Will be interpolated online via FlexGrid from native 2.0 x 2.5 resolution\n"
elif [[ ${met} = "geosit" ]]; then
    printf "  1. 4.0  x 5.0\n"
    printf "  2. 2.0  x 2.5\n"
    if [[ "${sim_name}" = "TransportTracers" ]]; then
	printf "  3. 0.5  x 0.625\n"
    fi
else
    printf "  1. 4.0  x 5.0\n"
    printf "  2. 2.0  x 2.5\n"
    printf "  3. 0.5  x 0.625\n"
    if [[ "x${met}" == "xgeosfp" ]]; then
	printf "  4. 0.25 x 0.3125\n"
    fi
fi

valid_res=0
while [ "${valid_res}" -eq 0 ]; do
    read -p "${USER_PROMPT}" res_num
    valid_res=1
    if [[ "x${res_num}" == "x1" ]]; then
	grid_res='4x5'
	RUNDIR_VARS+="$(cat ${gcdir}/run/shared/settings/4x5.txt)\n"
    elif [[ "x${res_num}" == "x2" ]]; then
	grid_res='2x25'
	RUNDIR_VARS+="$(cat ${gcdir}/run/shared/settings/2x25.txt)\n"
    elif [[ "x${res_num}" == "x3" ]]; then
	grid_res='05x0625'
	RUNDIR_VARS+="$(cat ${gcdir}/run/shared/settings/05x0625.txt)\n"
    elif [[ "x${res_num}" == "x4" ]]; then
	# Error check: Don't allow a 0.25 x 0.3125 MERRA-2 rundir.
	#  -- Melissa Sulprizio, Bob Yantosca (12 Sep 2023)
	if [[ "x${met}" == "xmerra2" ]]; then
	    valid_res=0
	    printf "Cannot create a MERRA-2 rundir at 0.25 x 0.3125 "
	    printf "resolution!\nPlease make another selection.\n"
	fi
	grid_res='025x03125'
	RUNDIR_VARS+="$(cat ${gcdir}/run/shared/settings/025x03125.txt)\n"
    else
	valid_res=0
	printf "Invalid horizontal resolution option.\n"
	printf "Please make another selection.\n"
    fi
done

if [[ "${met}" != "geosit" ]] && [[ "${grid_res}" = "05x0625" || "${grid_res}" = "025x03125" ]]; then
    printf "${thinline}Choose horizontal grid domain:${thinline}"
    printf "  1. Global\n"
    printf "  2. Asia\n"
    printf "  3. Europe\n"
    printf "  4. North America\n"
    printf "  5. Custom\n"

    valid_domain=0
    while [ "${valid_domain}" -eq 0 ]; do
	read -p "${USER_PROMPT}" domain_num
	valid_domain=1
	if [[ ${domain_num} = "1" ]]; then
	    RUNDIR_VARS+="$(cat ${gcdir}/run/shared/settings/global_grid.txt)\n"
	    RUNDIR_VARS+="RUNDIR_GRID_HALF_POLAR='true '\n"
	else
	    RUNDIR_VARS+="$(cat ${gcdir}/run/shared/settings/nested_grid.txt)\n"
	    if [[ ${domain_num} = "2" ]]; then
		RUNDIR_VARS+="RUNDIR_GRID_DOMAIN_NAME='AS'\n"
		grid_nest="AS"
	        if [[ ${grid_res} = "05x0625" ]]; then
	            RUNDIR_VARS+="RUNDIR_GRID_LON_RANGE='[ 60.0, 150.0]'\n"
		    RUNDIR_VARS+="RUNDIR_GRID_LAT_RANGE='[-11.0,  55.0]'\n"
		elif [[ ${grid_res} = "025x03125" ]]; then
	            RUNDIR_VARS+="RUNDIR_GRID_LON_RANGE='[ 70.0, 140.0]'\n"
		    RUNDIR_VARS+="RUNDIR_GRID_LAT_RANGE='[ 15.0,  55.0]'\n"
		fi
	    elif [[ ${domain_num} = "3" ]]; then
		RUNDIR_VARS+="RUNDIR_GRID_DOMAIN_NAME='EU'\n"
	        grid_nest="EU"
	        if [[ ${grid_res} = "05x0625" ]]; then
	            RUNDIR_VARS+="RUNDIR_GRID_LON_RANGE='[-30.0, 50.0]'\n"
		    RUNDIR_VARS+="RUNDIR_GRID_LAT_RANGE='[ 30.0, 70.0]'\n"
		elif [[ ${grid_res} = "025x03125" ]]; then
	            RUNDIR_VARS+="RUNDIR_GRID_LON_RANGE='[-15.0,  40.0 ]'\n"
		    RUNDIR_VARS+="RUNDIR_GRID_LAT_RANGE='[ 32.75, 61.25]'\n"
		fi
	    elif [[ ${domain_num} = "4" ]]; then
		RUNDIR_VARS+="RUNDIR_GRID_DOMAIN_NAME='NA'\n"
		grid_nest+="NA"
	        if [[ ${grid_res} = "05x0625" ]]; then
	            RUNDIR_VARS+="RUNDIR_GRID_LON_RANGE='[-140.0, -40.0]'\n"
		    RUNDIR_VARS+="RUNDIR_GRID_LAT_RANGE='[  10.0,  70.0]'\n"
		elif [[ ${grid_res} = "025x03125" ]]; then
	            RUNDIR_VARS+="RUNDIR_GRID_LON_RANGE='[-130.0,  -60.0]'\n"
		    RUNDIR_VARS+="RUNDIR_GRID_LAT_RANGE='[   9.75,  60.0]'\n"
		fi
	    elif [[ ${domain_num} = "5" ]]; then
		grid_nest="CU"
		RUNDIR_VARS+="RUNDIR_GRID_DOMAIN_NAME='custom'\n"
	        RUNDIR_VARS+="RUNDIR_GRID_LON_RANGE='MinLon MaxLon'\n"
	        RUNDIR_VARS+="RUNDIR_GRID_LAT_RANGE='MinLat MaxLat'\n"
	        printf "\n  -- You will need to manually set longitude and latitude"
		printf "\n     bounds in the Grid Menu of geoschem_config.yml!\n"
	    else
  		valid_domain=0
		printf "Invalid horizontal grid domain option. Try again.\n"
	    fi
        fi
    done
else
    RUNDIR_VARS+="$(cat ${gcdir}/run/shared/settings/global_grid.txt)\n"
    if [[ ${met} = "ModelE2.1" ]] || [[ ${met} = "ModelE2.2" ]]; then
        if [[ "$grid_res" == "4x5" ]]; then
	    RUNDIR_VARS+="RUNDIR_GRID_HALF_POLAR='true '\n"
	else
	    RUNDIR_VARS+="RUNDIR_GRID_HALF_POLAR='false'\n"
	fi
    else
	RUNDIR_VARS+="RUNDIR_GRID_HALF_POLAR='true '\n"
    fi
fi


RUNDIR_VARS+="$(cat ${shared_met_settings})\n"   # shared_met_settings needs to be included after RUNDIR_GRID_DIR is defined

# Set timesteps according to grid resolution
if [[ ${grid_res} = "05x0625" ]] || [[ ${grid_res} = "025x03125" ]]; then
    RUNDIR_VARS+="RUNDIR_TRANSPORT_TS='300'\n"
    RUNDIR_VARS+="RUNDIR_CHEMISTRY_TS='600'\n"
else
    if [[ ${sim_extra_option} =~ "TOMAS" ]]; then
	RUNDIR_VARS+="RUNDIR_TRANSPORT_TS='1800'\n"
	RUNDIR_VARS+="RUNDIR_CHEMISTRY_TS='3600'\n"
    else
	RUNDIR_VARS+="RUNDIR_TRANSPORT_TS='600'\n"
	RUNDIR_VARS+="RUNDIR_CHEMISTRY_TS='1200'\n"
    fi
fi

#-----------------------------------------------------------------
# Is International Date Line an edge or midpoint?
#-----------------------------------------------------------------

if [[ ${met} = "ModelE2.1" ]] || [[ ${met} = "ModelE2.2" ]] ; then
    if [[ "$grid_res" == "2x25" ]]; then
	# Native GISS fine resolution
	RUNDIR_VARS+="RUNDIR_CENTER_LON_180='false'\n"
    else
        # FlexGrid re-gridded resolutions
	RUNDIR_VARS+="RUNDIR_CENTER_LON_180='true '\n"
    fi
else
    # All GMAO products
    RUNDIR_VARS+="RUNDIR_CENTER_LON_180='true '\n"
fi

#----------------------------------------------------------------
# Horizontal resolution-dependent settings
#-----------------------------------------------------------------

if [[ ${met} = "ModelE2.1" ]]; then
    if [[ "$runid" == "E213f10aF40oQ40nudge" ]]; then
        if [[ "$grid_res" ==  "4x5" ]]; then
	    RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='0.00474046'\n"
        elif [[ "$grid_res" == "2x25" ]]; then
	    RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='0.00243979'\n"
        elif [[ "$grid_res" == "05x0625" ]]; then
	    RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='0.00276896'\n"
        elif [[ "$grid_res" == "025x03125" ]]; then
	    RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='0.00254319'\n"
	else
	    RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='-999.0e0'\n"
  	fi
    else
        if [[ "$grid_res" ==  "4x5" ]]; then
            RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='0.03564873'\n"
        elif [[ "$grid_res" == "2x25" ]]; then
	    RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='0.01050036'\n"
        elif [[ "$grid_res" == "05x0625" ]]; then
	    RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='0.01340854'\n"
        elif [[ "$grid_res" == "025x03125" ]]; then
	    RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='0.01066495'\n"
	else
	    RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='-999.0e0'\n"
	fi
    fi
else
    RUNDIR_VARS+="RUNDIR_GISS_RES='not_used'\n"
    # Use GEOS-FP values as placeholders for GEOS-IT until parameters derived
    if [[ "x${sim_name}" == "xfullchem" || "x${sim_name}" == "xaerosol" ]]; then
	if [[ "x${met}" == "xgeosfp" && "x${grid_res}" == "x4x5" ]]; then
	    RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='8.3286e-4'\n"
	elif [[ "x${met}" == "xgeosfp" && "x${grid_res}" == "x2x25" ]]; then
	    RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='5.0416e-4'\n"
	elif [[ "x${met}" == "xmerra2" && "x${grid_res}" == "x4x5" ]]; then
	    RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='7.8533e-4'\n"
	elif [[ "x${met}" == "xmerra2" && "x${grid_res}" == "x2x25" ]]; then
	    RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='4.7586e-4'\n"
	elif [[ "x${met}" == "xgeosit" && "x${grid_res}" == "x4x5" ]]; then
	    RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='8.3286e-4'\n"
	elif [[ "x${met}" == "xgeosit" && "x${grid_res}" == "x2x25" ]]; then
	    RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='5.0416e-4'\n"
	else
	    RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='-999.0e0'\n"
	fi
    else
	RUNDIR_VARS+="RUNDIR_DUSTDEAD_TF='-999.0e0'\n"
    fi
fi

#-----------------------------------------------------------------
# Ask user to select vertical resolution
#-----------------------------------------------------------------
printf "${thinline}Choose number of levels:${thinline}"

if [[ ${met} = "geosfp" ]] || [[ ${met} = "merra2" || ${met} = "geosit" ]]; then
    printf "  1. 72 (native)\n"
    printf "  2. 47 (reduced)\n"

    valid_lev=0
    while [ "${valid_lev}" -eq 0 ]; do
        read -p "${USER_PROMPT}" lev_num
        valid_lev=1
        if [[ ${lev_num} = "1" ]]; then
            RUNDIR_VARS+="RUNDIR_GRID_NLEV='72'\n"
        elif [[ ${lev_num} = "2" ]]; then
	    grid_lev="47L"
	    RUNDIR_VARS+="RUNDIR_GRID_NLEV='47'\n"
        else
            valid_lev=0
            printf "Invalid vertical resolution option. Try again.\n"
        fi
    done
fi

if [[ ${met} = "ModelE2.1" ]]; then
    printf "  1. 40 (native)\n"
    valid_lev=0
    while [ "${valid_lev}" -eq 0 ]; do
        read -p "${USER_PROMPT}" lev_num
        valid_lev=1
        if [[ ${lev_num} = "1" ]]; then
            RUNDIR_VARS+="RUNDIR_GRID_NLEV='40'\n"
        else
            valid_lev=0
            printf "Invalid vertical resolution option. Try again.\n"
        fi
    done
fi

if [[ ${met} = "ModelE2.2" ]]; then
    printf "  1. 102 (native)\n"
    printf "  2. 74 (reduced)\n"
    valid_lev=0
    while [ "${valid_lev}" -eq 0 ]; do
        read -p "${USER_PROMPT}" lev_num
        valid_lev=1
        if [[ ${lev_num} = "1" ]]; then
            RUNDIR_VARS+="RUNDIR_GRID_NLEV='102'\n"
        elif [[ ${lev_num} = "2" ]]; then
            RUNDIR_VARS+="RUNDIR_GRID_NLEV='74'\n"
        else
            valid_lev=0
            printf "Invalid vertical resolution option. Try again.\n"
        fi
    done
fi

#-----------------------------------------------------------------
# Ask user to define path where directory will be created
#-----------------------------------------------------------------
printf "${thinline}Enter path where the run directory will be created:${thinline}"
valid_path=0
while [ "$valid_path" -eq 0 ]; do
    read -e -p "${USER_PROMPT}" rundir_path

    # Test for quitting
    if [[ "x${rundir_path}" == "xq" ]]; then
	printf "\nExiting.\n"
	exit 1
    fi

    # Expand $rundir_path to an absolute path.
    # Also replace ~ with $HOME.
    rundir_path="${rundir_path/#\~/$HOME}"
    rundir_path=$(realpath "${rundir_path}")
    echo "Expanding to ${rundir_path}"

    # If this is just a new directory within an existing one,
    # give the user the option to proceed
    if [[ ! -d ${rundir_path} ]]; then
        if [[ -d $(dirname ${rundir_path} ) ]]; then
            printf "\nWarning: ${rundir_path} does not exist,\nbut the parent directory does.\nWould you like to make this directory? (y/n/q)\n"
            read -p "${USER_PROMPT}" mk_rundir
            if [[ "x${mk_rundir}" == "xy" ]]; then
                mkdir $rundir_path
	    elif [[ "x${mk_rundir}" == "xq" ]]; then
		printf "\nExiting.\n"
		exit 1
            fi
        fi
    fi

    # Ask user to supply a new path again
    if [[ ! -d ${rundir_path} ]]; then
        printf "\nERROR: ${rundir_path} does not exist. Enter a new path or hit q to quit.\n"
    else
	valid_path=1
    fi
done

#-----------------------------------------------------------------
# Ask user to define run directory name if not passed as argument
#-----------------------------------------------------------------
if [ -z "$1" ]; then
    printf "${thinline}Enter run directory name, or press return to use default:\n\n"
    printf "NOTE: This will be a subfolder of the path you entered above.${thinline}"
    read -e -p "${USER_PROMPT}" rundir_name
    if [[ -z "${rundir_name}" ]]; then
	grid_display="${grid_res}"
        [[ "x${grid_nest}" != "x" ]] && grid_display+="_${grid_nest}"
	[[ "x${grid_lev}"  != "x" ]] && grid_display+="_${grid_lev}"
	if [[ "x${sim_extra_option}" = "xnone" ]]; then
	    rundir_name=gc_${grid_display}_${met}_${sim_name}
	else
	    rundir_name=gc_${grid_display}_${met}_${sim_name}_${sim_extra_option}
	fi
	printf "  -- Using default directory name ${rundir_name}"
    fi
else
    rundir_name=$1
fi

#-----------------------------------------------------------------
# Ask user for a new run directory name if specified one exists
#-----------------------------------------------------------------
rundir=${rundir_path}/${rundir_name}
valid_rundir=0
while [ "${valid_rundir}" -eq 0 ]; do
    if [[ -d ${rundir} ]]; then
	printf "\nWARNING: ${rundir} already exists.\n"
        printf "Enter a different run directory name, or q to quit:\n"
	read -e -p "${USER_PROMPT}" new_rundir
	if [[ ${new_rundir} = "q" ]]; then
	    printf "Exiting.\n"
	    exit 1
	else
	    rundir=${rundir_path}/${new_rundir}
	fi
    else
        valid_rundir=1
    fi
done

#-----------------------------------------------------------------
# Create run directory
#-----------------------------------------------------------------

mkdir -p ${rundir}
mkdir -p ${rundir}/Restarts

# Copy run directory files and subdirectories
cp ${gcdir}/run/shared/cleanRunDir.sh       ${rundir}
cp ${gcdir}/run/shared/download_data.py     ${rundir}
cp ${gcdir}/run/shared/download_data.yml    ${rundir}
cp ./getRunInfo                             ${rundir}
cp ./archiveRun.sh                          ${rundir}
cp ./README.md                              ${rundir}
cp ./gitignore                              ${rundir}/.gitignore

# Use data downloader that points to GCAP2 restart files
if [[ ${met} = "ModelE2.1" ]] || [[ ${met} = "ModelE2.2" ]]; then
  cp ${gcdir}/run/shared/download_data.gcap2.40L.yml ${rundir}/download_data.yml
fi

# Copy the OH metrics Python script to the rundir (fullchem/CH4 only)
if [[ "x${sim_name}" == "xfullchem" || "x${sim_name}" == "xCH4" ]]; then
    cp -r ${gcdir}/run/shared/metrics.py  ${rundir}
    chmod 744 ${rundir}/metrics.py
fi

# Set permissions
chmod 744 ${rundir}/cleanRunDir.sh
chmod 744 ${rundir}/archiveRun.sh

# Copy species database.  NOTE: For the new Hg simulation via KPP, we need
# to copy species_database_hg.yml to the rundir and rename it to
# species_database.yml.  This is because the Hg simulation has several
# inactive species that are active in the other simulations, and this
# causes a conflict.  Work out a better solution later.
#  -- Bob Yantosca, 10 Dec 2021
if [[ "x${sim_name}" == "xHg" ]]; then
    cp -r ${gcdir}/run/shared/species_database_hg.yml ${rundir}/species_database.yml
else
    cp -r ${gcdir}/run/shared/species_database.yml ${rundir}
fi

# Append APM or TOMAS species to species database
# Also copy APM input files to the run directory
if [[ ${sim_extra_option} =~ "TOMAS" ]]; then
    cat ${gcdir}/run/shared/species_database_tomas.yml >> ${rundir}/species_database.yml
elif [[ ${sim_extra_option} =~ "APM" ]]; then
    cat ${gcdir}/run/shared/species_database_apm.yml >> ${rundir}/species_database.yml
    cp ${gcdir}/run/shared/apm_tmp.dat ${rundir}/apm_tmp.dat
    cp ${gcdir}/run/shared/input.apm   ${rundir}/input.apm
fi

# If benchmark simulation, put run script in directory
if [[ "x${sim_extra_option}" == "xbenchmark" ]]; then
    scriptDir="./runScriptSamples/operational_examples/harvard_cannon"
    cp ${scriptDir}/geoschem.benchmark.run ${rundir}
    chmod 744 ${rundir}/geoschem.benchmark.run
fi

# Create symbolic link to code directory
ln -s ${wrapperdir} ${rundir}/CodeDir
ln -s ${wrapperdir}/run/GCClassic/runScriptSamples ${rundir}/runScriptSamples

#--------------------------------------------------------------------
# Navigate to run directory and set up input files
#--------------------------------------------------------------------
cd ${rundir}

# Special handling for start/end date based on simulation so that
# start year/month/day matches default initial restart file.
if [[ "x${sim_name}" == "xHg"     ||
      "x${sim_name}" == "xCH4"    ||
      "x${sim_name}" == "xCO2"    ||
      "x${sim_name}" == "xtagCO"  ||
      "x${sim_name}" == "xcarbon" ||
      "x${sim_name}" == "xTransportTracers" ]]; then
    startdate='20190101'
    enddate='20190201'
elif [[ "x${sim_name}" == "xmetals" ]]; then
    startdate='20110101'
    enddate='20110201'
else
    startdate='20190701'
    enddate='20190801'
fi

if [[ ${met} = "ModelE2.1" ]] || [[ ${met} = "ModelE2.2" ]]; then
    if [[ "$scenario" == "HIST" ]]; then
	startdate='20050701'
	enddate='20050801'
    else
	startdate='20900701'
	enddate='20900801'
    fi
fi
RUNDIR_VARS+="RUNDIR_SIM_START_DATE=$startdate\n"
RUNDIR_VARS+="RUNDIR_SIM_END_DATE=$enddate\n"
RUNDIR_VARS+="RUNDIR_SIM_START_TIME='000000'\n"
RUNDIR_VARS+="RUNDIR_SIM_END_TIME='000000'\n"

# Use monthly diagnostics by default
RUNDIR_VARS+="RUNDIR_HIST_TIME_AVG_DUR='00000100 000000'\n"
RUNDIR_VARS+="RUNDIR_HIST_TIME_AVG_FREQ='00000100 000000'\n"
RUNDIR_VARS+="RUNDIR_HIST_INST_DUR='00000100 000000'\n"
RUNDIR_VARS+="RUNDIR_HIST_INST_FREQ='00000100 000000'\n"
RUNDIR_VARS+="RUNDIR_HIST_MONTHLY_DIAG='1'\n"

# Turn on GEOS-Chem timers for benchmark simulations
if [[ "${sim_extra_option}" == "benchmark" ]]; then
    RUNDIR_VARS+="RUNDIR_USE_GCCLASSIC_TIMERS='true '\n"
else
    RUNDIR_VARS+="RUNDIR_USE_GCCLASSIC_TIMERS='false'\n"
fi

# Assign appropriate file paths and settings in HEMCO_Config.rc
if [[ ${met} = "ModelE2.1" ]]; then
    RUNDIR_VARS+="RUNDIR_DUSTDEAD_EXT='on '\n"
    RUNDIR_VARS+="RUNDIR_MEGAN_EXT='on '\n"
    RUNDIR_VARS+="RUNDIR_SEASALT_EXT='on '\n"
    RUNDIR_VARS+="RUNDIR_SOILNOX_EXT='on '\n"
    RUNDIR_VARS+="RUNDIR_OFFLINE_DUST='false'\n"
    RUNDIR_VARS+="RUNDIR_OFFLINE_BIOVOC='false'\n"
    RUNDIR_VARS+="RUNDIR_OFFLINE_SEASALT='false'\n"
    RUNDIR_VARS+="RUNDIR_OFFLINE_SOILNOX='false'\n"
    RUNDIR_VARS+="RUNDIR_TOMAS_SEASALT='off'\n"
    RUNDIR_VARS+="RUNDIR_TOMAS_DUSTDEAD='off'\n"
    RUNDIR_VARS+="$(cat ${gcdir}/run/shared/settings/gcap2_hemco.txt)\n"
else
    if [[ "${sim_extra_option}" == "benchmark" ]]; then
	RUNDIR_VARS+="RUNDIR_DUSTDEAD_EXT='on '\n"
	RUNDIR_VARS+="RUNDIR_MEGAN_EXT='on '\n"
	RUNDIR_VARS+="RUNDIR_SEASALT_EXT='on '\n"
	RUNDIR_VARS+="RUNDIR_SOILNOX_EXT='on '\n"
	RUNDIR_VARS+="RUNDIR_OFFLINE_DUST='false'\n"
	RUNDIR_VARS+="RUNDIR_OFFLINE_BIOVOC='false'\n"
	RUNDIR_VARS+="RUNDIR_OFFLINE_SEASALT='false'\n"
	RUNDIR_VARS+="RUNDIR_OFFLINE_SOILNOX='false'\n"
	RUNDIR_VARS+="RUNDIR_TOMAS_SEASALT='off'\n"
	RUNDIR_VARS+="RUNDIR_TOMAS_DUSTDEAD='off'\n"
    else
	if [[ "${sim_extra_option}" == "marinePOA" ]]; then
	    RUNDIR_VARS+="RUNDIR_SEASALT_EXT='on '\n"
	    RUNDIR_VARS+="RUNDIR_OFFLINE_SEASALT='false'\n"
	    RUNDIR_VARS+="RUNDIR_TOMAS_SEASALT='off'\n"
	else
	    RUNDIR_VARS+="RUNDIR_SEASALT_EXT='off'\n"
	    if [[ ${sim_extra_option} =~ "TOMAS" ]]; then
		RUNDIR_VARS+="RUNDIR_TOMAS_SEASALT='on '\n"
		RUNDIR_VARS+="RUNDIR_OFFLINE_SEASALT='false'\n"
	    else
		RUNDIR_VARS+="RUNDIR_TOMAS_SEASALT='off'\n"
		RUNDIR_VARS+="RUNDIR_OFFLINE_SEASALT='true '\n"
	    fi
	fi
	if [[ ${sim_extra_option} =~ "TOMAS" ]]; then
	    RUNDIR_VARS+="RUNDIR_TOMAS_DUSTDEAD='on '\n"
	    RUNDIR_VARS+="RUNDIR_OFFLINE_DUST='false'\n"
	else
	    RUNDIR_VARS+="RUNDIR_TOMAS_DUSTDEAD='off'\n"
	    RUNDIR_VARS+="RUNDIR_OFFLINE_DUST='true '\n"
	fi
	RUNDIR_VARS+="RUNDIR_DUSTDEAD_EXT='off'\n"
	RUNDIR_VARS+="RUNDIR_MEGAN_EXT='off'\n"
	RUNDIR_VARS+="RUNDIR_SOILNOX_EXT='off'\n"
	RUNDIR_VARS+="RUNDIR_OFFLINE_BIOVOC='true '\n"
	RUNDIR_VARS+="RUNDIR_OFFLINE_SOILNOX='true '\n"
    fi
    RUNDIR_VARS+="$(cat ${gcdir}/run/shared/settings/gmao_hemco.txt)\n"
fi

#--------------------------------------------------------------------
# Replace settings in config files with RUNDIR variables
#--------------------------------------------------------------------

# Define a subdirectory for rundir configuration files
rundir_config_dirname=CreateRunDirLogs
mkdir -p ${rundir}/${rundir_config_dirname}

# Save RUNDIR variables to a file in the rundirConfig folder
rundir_config_logname=rundir_vars.txt
rundir_config_log=${rundir}/${rundir_config_dirname}/${rundir_config_logname}
echo -e "$RUNDIR_VARS" > ${rundir_config_log}
#sort -o ${rundir_config_log} ${rundir_config_log}

# Initialize run directory
# NOTE: This also copies configuration files to the run directory!
${srcrundir}/init_rd.sh ${rundir_config_log}

#--------------------------------------------------------------------
# Print run directory setup info to screen
#--------------------------------------------------------------------

printf "\n  -- See ${rundir_config_dirname}/${rundir_config_logname} for summary of default run directory settings"
printf "\n  -- This run directory has been set up to start on ${startdate}"
printf "\n  -- A restart file for this date has been copied to the Restarts subdirectory"
printf "\n  -- You may add more restart files using format GEOSChem.Restart.YYYYMMDD_HHmmz.nc4"
printf "\n  -- Change simulation start and end dates in configuration file geoschem_config.yml"
printf "\n  -- Default frequency and duration of diagnostics are set to monthly"
printf "\n  -- Modify diagnostic settings in HISTORY.rc and HEMCO_Config.rc\n"

if [[ "x${nested_sim}" == "xT" ]]; then
    printf "\n  -- Nested-grid simulations use global high-reoslution met fields"
    printf "\n     by default. To improve run time, you may choose to use cropped"
    printf "\n     met fields by modifying the file paths and names in HEMCO_Config.rc"
    printf "\n     to include the region string (e.g. 'AS', 'EU', 'NA').\n"
fi

#--------------------------------------------------------------------
# Copy sample restart file to run directory
# Bash functions used here are from ./setupForRestarts.sh
#--------------------------------------------------------------------

# Parse the download_data.yml file, which returns variable declarations
# prefixed by "RUNDIR_", such as: RUNDIR_restarts_root="GEOSCHEM_RESTARTS"
if [[ "x${met}" == "xModelE2.1" || "x${met}" == "xModelE2.2" ]]; then
    config_file="${gcdir}/run/shared/download_data.gcap2.40L.yml"
else
    config_file="${gcdir}/run/shared/download_data.yml"
fi
config=$(parseYaml "${config_file}" "RUNDIR_")

# Export environment variables to be used by rundir scripts
for var in ${config[@]}; do
    setRestartEnvVar "${var}"
done

# Root paths for restarts
# Check the Linux Kernel version to see if we are on the AWS cloud.
# If we are, define the command to copy the restart file from s3://gcgrid
is_aws=$(uname -r | grep aws)
rst_root=$(getRemoteRoot "${is_aws}")
s3_cp=$(getS3CopyCmd "${is_aws}")
loc_root="${rundir}/Restarts"

# Copy the proper restart file to the run directory Restarts/ folder
copyRestartToRunDir "${sim_name}" "${sim_extra_option}" \
		    "${rst_root}" "${loc_root}"

# Change time cycle flags in HEMCO_Config.rc for those simulations
# in which the restart files do not contain all species
setEFYOtoEYinHemcoConfig "${sim_name}" "${sim_extra_option}"

# Unset environment variables used by rundir scripts
for var in ${config[@]}; do
    unsetRestartEnvVar "${var}"
done

#--------------------------------------------------------------------
# Other setup tasks
#--------------------------------------------------------------------

# Call function to setup configuration files with settings common between
# GEOS-Chem Classic and GCHP. This script mainly now adds species to
# input_options.yml and modifies diagnostic output based on simulation type.
if [[ "x${sim_name}" = "xfullchem" ]]; then
    set_common_settings "${sim_extra_option}" "GCClassic"
fi 

# If necessary, edit config files for a carbon single species simulation
if [[ "x${sim_name}" == "xcarbon" ]]; then
    if [[ "x${sim_extra_option}" != "xnone" ]]; then
	singleCarbonSpecies "${sim_extra_option}" "${rundir}"
    fi
fi

#--------------------------------------------------------------------
# Navigate back to source code directory
#--------------------------------------------------------------------
cd ${srcrundir}

#----------------------------------------------------------------------
# Archive repository version in run directory file rundir.version
#----------------------------------------------------------------------

# Get info about the current commit in geoschem/geos-chem
cd ${gcdir}
remote_url=$(git config --get remote.origin.url)
code_branch=$(git rev-parse --abbrev-ref HEAD)
last_commit=$(git log -n 1 --pretty=format:"%s")
commit_date=$(git log -n 1 --pretty=format:"%cd")
commit_user=$(git log -n 1 --pretty=format:"%cn")
commit_hash=$(git log -n 1 --pretty=format:"%h")
cd ${srcrundir}

# Write commit info to a version log
version_log=${rundir}/${rundir_config_dirname}/rundir.version
printf   " This run directory was created with:"        >  ${version_log}
printf "\n ${srcrundir}/createRunDir.sh.\n"             >> ${version_log}
printf "\n GEOS-Chem repository version information:\n" >> ${version_log}
printf "\n  Remote URL: ${remote_url}"                  >> ${version_log}
printf "\n  Branch: ${code_branch}"                     >> ${version_log}
printf "\n  Commit: ${last_commit}"                     >> ${version_log}
printf "\n  Date: ${commit_date}"                       >> ${version_log}
printf "\n  User: ${commit_user}"                       >> ${version_log}
printf "\n  Hash: ${commit_hash}\n"                     >> ${version_log}

#-----------------------------------------------------------------
# Ask user whether to track run directory changes with git
#-----------------------------------------------------------------
printf "${thinline}Do you want to track run directory changes with git? (y/n)${thinline}"
valid_response=0
while [ "$valid_response" -eq 0 ]; do
    read -p "${USER_PROMPT}" enable_git
    if [[ "x${enable_git}" == "xy" ]]; then
	cd ${rundir}
	printf "\n\nChanges to the following run directory files are tracked by git:\n\n" >> ${version_log}
	printf "\n"
	git init
	git add *.rc *.sh *.yml *.py geoschem_config.yml getRunInfo
	[[ -f geoschem.benchmark.run         ]] && git add geoschem.benchmark.run
	[[ -f geoschem.run                   ]] && git add geoschem.run
	[[ -f HEMCO_Config.rc.gmao_metfields ]] && git add HEMCO_Config.rc.gmao_metfields
	[[ -f ${rundir_config_log}           ]] && git add ${rundir_config_log}
	printf " " >> ${version_log}
	git commit -m "Initial run directory" >> ${version_log}
	cd ${srcrundir}
	valid_response=1
    elif [[ ${enable_git} = "n" ]]; then
	valid_response=1
    else
	printf "Input not recognized. Try again.\n"
    fi
done

#-----------------------------------------------------------------
# Ask user whether to build the KPP-standalone box model
#-----------------------------------------------------------------
enable_kppsa=""
if [[ "x${sim_name}" == "xfullchem" ]]; then
    printf "${thinline}Do you want to build the KPP-Standalone Box Model? (y/n)${thinline}"
    valid_response=0
    while [ "$valid_response" -eq 0 ]; do
	read -p "${USER_PROMPT}" enable_kppsa
	if [[ "x${enable_kppsa}" == "xy" ]]; then
	    cp -r ${gcdir}/run/shared/kpp_standalone_interface.yml  ${rundir}
	    chmod 644 ${rundir}/kpp_standalone_interface.yml
	    valid_response=1
	elif [[ "x${enable_kppsa}" = "xn" ]]; then
	    valid_response=1
	else
	    printf "Input not recognized. Try again.\n"
	fi
    done
fi

#---------------------------------------------------------------------------
# Add reminders to compile with CMake options for simulations that need them
#---------------------------------------------------------------------------
hdr="\n>>>> REMINDER: You must compile with options:"
ftr="<<<<\n"
EXTRA_CMAKE_OPTIONS=""
[[ "x${sim_name}" == "xcarbon" ]] && EXTRA_CMAKE_OPTIONS="-DMECH=carbon "
[[ "x${sim_name}" == "xHg"     ]] && EXTRA_CMAKE_OPTIONS="-DMECH=Hg -DFASTJX=y "
if [[ "x${sim_name}" == "xfullchem" ]]; then
    [[ "x${sim_extra_option}" == "xAPM"     ]] && EXTRA_CMAKE_OPTIONS="-DAPM=y "
    [[ "x${sim_extra_option}" == "xRRTMG"   ]] && EXTRA_CMAKE_OPTIONS="-DRRTMG=y "
    [[ "x${sim_extra_option}" == "xTOMAS15" ]] && EXTRA_CMAKE_OPTIONS="-DTOMAS=y -DTOMAS_BINS=15 "
    [[ "x${sim_extra_option}" == "xTOMAS40" ]] && EXTRA_CMAKE_OPTIONS="-DTOMAS=y -DTOMAS_BINS=40 "
fi
[[ "x${enable_kppsa}" == "xy" ]] && EXTRA_CMAKE_OPTIONS+="-DKPPSA=y"

# Add to RUNDIR_VARS
RUNDIR_VARS+="EXTRA_CMAKE_OPTIONS=${EXTRA_CMAKE_OPTIONS}"

# Print a reminder to compile with extra CMake options, if necessary
[[ "x${EXTRA_CMAKE_OPTIONS}" != "x" ]] && printf "${hdr} ${EXTRA_CMAKE_OPTIONS} ${ftr}"

#---------------------------------------------------------------------------
# Create build directory README file
#---------------------------------------------------------------------------
mkdir -p "${rundir}/build"
msg="To build GEOS-Chem, type:\n\n"
msg+="$ cmake ../CodeDir\n"
msg+="$ cmake . -DRUNDIR=.. ${EXTRA_CMAKE_OPTIONS}\n"
msg+="$ make -j\n"
msg+="$ make install\n"
printf "${msg}" > ${rundir}/build/README
unset msg

#-----------------------------------------------------------------
# Add the version info to the top of the rundir configuration log
#-----------------------------------------------------------------

# Add a caveat that these rundir settings only go with this commit
printf "\n\n IMPORTANT: ONLY USE THESE RUNDIR SETTINGS WITH THIS COMMIT!\n" >> ${version_log}

# Add a "# " characters to the front of each line so we can use
# it as a comment heading for ${rundir_config_logname}
sed 's/^/# /' ${version_log} > tmp.txt
mv tmp.txt ${version_log}

# Add the version log to the top of the rundir config log
cat ${version_log} ${rundir_config_log} > tmp.txt
mv tmp.txt ${rundir_config_log}

# Remove the version log
rm -rf ${version_log}

# Save the updated rundir_vars file to the git repo
if [[ "x${enable_git}" == "xy" ]]; then
    if [[ -f ${rundir_config_log} ]]; then
	cd ${rundir}
	git add ${rundir_config_log}
	git commit -m "Update header of ${rundir_config_dirname}/${rundir_config_logname}" > /dev/null
	cd ${srcrundir}
    fi
fi

#-----------------------------------------------------------------
# Create and populate warnings file
#-----------------------------------------------------------------
if [[ $met == "geosfp" ]]; then
   echo -e ${fp_msg} > ${rundir}/warnings.txt
fi

#-----------------------------------------------------------------
# Done!
#-----------------------------------------------------------------
printf "\nCreated ${rundir}\n"

exit 0
