#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: createRunDir.sh
#
# !DESCRIPTION: Creates a GCHP run directory.
#\\
#\\
# !CALLING SEQUENCE:
#  ./createRunDir.sh [rundirname]
#
# !REMARKS:
#  If optional run directory name argument is not passed then the user
#  will be prompted to enter a name interactively, or choose to use the
#  default name gchp_{met}_{sim_name}_{sim_extra_option}.
#
# !REVISION HISTORY:
#  Initial version: E. Lundgren,10/5/2018
#  See the subsequent Git history with the gitk browser!
#------------------------------------------------------------------------------
#BOC

# Directory w/ GCHP rundir scripts (i.e. this directory)
srcrundir=$(pwd -P)
cd ${srcrundir}

# GEOS-Chem "science codebase" directory
cd ../..
gcdir=$(pwd -P)

# GCHP "wrapper" directory
cd ../../../..
wrapperdir=$(pwd -P)

# Return to 
cd ${srcrundir}

# Source common bash functions from scripts in the run/shared folder
. ${gcdir}/run/shared/setupConfigFiles.sh      # Config file editing
. ${gcdir}/run/shared/newUserRegistration.sh   # 1st-time user registration
. ${gcdir}/run/shared/singleCarbonSpecies.sh   # Single carbon species setup

# Initialize run directory variables
RUNDIR_VARS=""
RUNDIR_VARS+="RUNDIR_GC_MODE='GCHP'\n"

# Define separator lines
thickline="\n===========================================================\n"
thinline="\n-----------------------------------------------------------\n"

printf "${thickline}GCHP RUN DIRECTORY CREATION${thickline}"

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
[[ -z "${GC_USER_REGISTERED}" ]] && registerNewUser "gchp"

#-----------------------------------------------------------------
# Ask user to select simulation type
#-----------------------------------------------------------------
printf "${thinline}Choose simulation type:${thinline}"
printf "   1. Full chemistry\n"
printf "   2. TransportTracers\n"
printf "   3. CO2 w/ CMS-Flux emissions\n"
printf "   4. Tagged O3\n"
printf "   5. Carbon\n"

valid_sim=0
while [ "${valid_sim}" -eq 0 ]; do
    read -p "${USER_PROMPT}" sim_num
    valid_sim=1
    if [[ ${sim_num} = "1" ]]; then
	sim_name=fullchem
    elif [[ ${sim_num} = "2" ]]; then
	sim_name=TransportTracers
    elif [[ ${sim_num} = "3" ]]; then
	sim_name=CO2
    elif [[ ${sim_num} = "4" ]]; then
	sim_name=tagO3
    elif [[ ${sim_num} = "5" ]]; then
	sim_name=carbon
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
if [[ ${sim_extra_option} == "benchmark"  ]] || \
   [[ ${sim_extra_option} =~ "complexSOA" ]] || \
   [[ ${sim_extra_option} == "APM"        ]]; then
    RUNDIR_VARS+="RUNDIR_COMPLEX_SOA='true'\n"
    if [[ ${sim_extra_option} == "complexSOA_SVPOA" ]]; then
	RUNDIR_VARS+="RUNDIR_SVPOA='true'\n"
    else
	RUNDIR_VARS+="RUNDIR_SVPOA='false'\n"
    fi
else
    RUNDIR_VARS+="RUNDIR_COMPLEX_SOA='false'\n"
    RUNDIR_VARS+="RUNDIR_SVPOA='false'\n"
fi

if [[ ${sim_extra_option} == "aciduptake" ]]; then
    RUNDIR_VARS+="RUNDIR_DUSTALK_EXT='on '\n"
    RUNDIR_VARS+="RUNDIR_ACID_UPTAKE='true'\n"
else
    RUNDIR_VARS+="RUNDIR_DUSTALK_EXT='off'\n"
    RUNDIR_VARS+="RUNDIR_ACID_UPTAKE='false'\n"
fi

if [[ ${sim_extra_option} == "marinePOA" ]]; then
    RUNDIR_VARS+="RUNDIR_MARINE_POA='true'\n"
else
    RUNDIR_VARS+="RUNDIR_MARINE_POA='false'\n"
fi

if [[ ${sim_extra_option} == "RRTMG" ]]; then
    RUNDIR_VARS+="RUNDIR_RRTMG_OPTS='true'\n"
    RUNDIR_VARS+="RUNDIR_USE_RRTMG='true '\n"
else
    RUNDIR_VARS+="RUNDIR_RRTMG_OPTS='false'\n"
    RUNDIR_VARS+="RUNDIR_USE_RRTMG='false'\n"
fi

if [[ ${sim_extra_option} =~ "TOMAS" ]]; then
    RUNDIR_VARS+="RUNDIR_USE_NLPBL='false'\n"
    RUNDIR_VARS+="RUNDIR_USE_ONLINE_O3='false'\n"
else
    RUNDIR_VARS+="RUNDIR_USE_NLPBL='true'\n"
    RUNDIR_VARS+="RUNDIR_USE_ONLINE_O3='true'\n"
fi

# NOTE: Fullchem benchmarks use the climatological volcano emissions!
if [[ "x${sim_name}" == "xfullchem" ]]; then
    RUNDIR_VARS+="RUNDIR_VOLC_CLIMATOLOGY='\$ROOT/VOLCANO/v2024-04/so2_volcanic_emissions_CARN_v202401.degassing_only.rc'\n"

    if [[ "x${sim_extra_option}" == "xbenchmark" ]]; then
	RUNDIR_VARS+="RUNDIR_VOLC_TABLE='\$ROOT/VOLCANO/v2024-04/so2_volcanic_emissions_CARN_v202401.degassing_only.rc'\n"
    else
	RUNDIR_VARS+="RUNDIR_VOLC_TABLE='\$ROOT/VOLCANO/v2024-04/\$YYYY/\$MM/so2_volcanic_emissions_Carns.\$YYYY\$MM\$DD.rc'\n"
    fi
fi

#-----------------------------------------------------------------
# Ask user to select meteorology source
#-----------------------------------------------------------------
printf "${thinline}Choose meteorology source:${thinline}"
printf "  1. MERRA-2 (Recommended)\n"
printf "  2. GEOS-FP\n"
printf "  3. GEOS-IT\n"

metSettingsDir=${gcdir}/run/shared/settings

valid_met=0
while [ "${valid_met}" -eq 0 ]; do
    read -p "${USER_PROMPT}" met_num
    valid_met=1

    if [[ ${met_num} = "1" ]]; then

	met="merra2"
	RUNDIR_VARS+="$(cat ${metSettingsDir}/merra2.txt)\n"

    elif [[ ${met_num} = "2" ]]; then

       	met="geosfp"

	# Print warning about GEOS-FP and require user to acknowledge it.
	fp_msg="WARNING: The convection scheme used to generate archived GEOS-FP meteorology \nfiles changed from RAS to Grell-Freitas starting June 1 2020 with impact on \nvertical transport. Discussion and analysis of the impact is available at \ngithub.com/geoschem/geos-chem/issues/1409. To fix this issue, different GEOS-Chem \nconvection schemes are called based on simulation start time. This ensures \ncomparability in GEOS-Chem runs using GEOS-FP fields generated using the RAS \nconvection scheme and fields generated using Grell-Freitas, but only if the \nsimulation does not cross the June 1 2020 boundary. We therefore recommend \nsplitting up GEOS-FP runs in time such that a single simulation does not span \nthis date. For example, configure one run to end on June 1 2020 and then use \nits output restart to start another run on June 1. Alternatively consider using \nMERRA2 which was entirely generated with RAS, or GEOS-IT which was entirely \ngenerated with Grell-Freitas. If you wish to use a GEOS-FP meteorology year \ndifferent from your simulation year please create a GEOS-Chem GitHub issue for \nassistance to avoid accidentally using zero convective precipitation flux.\n"
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

        # Ask user to specify processed or raw files
        printf "${thinline}Choose meteorology file type:${thinline}"
	printf "  1. 0.25x0.3125 processed files from the GEOS-Chem data archive\n"
	printf "  2. 0.25x0.3125 raw files that you downloaded from NASA GMAO\n"
	printf "\nNOTE: Processed files have unused variables removed, are vertically flipped,\n"
	printf "and are concatenated into daily files. Using them speeds up GCHP runs.\n\n"
	valid_response=0
	while [ "${valid_response}" -eq 0 ]; do
	    valid_response=1
	    read -p "${USER_PROMPT}" response
	    if [[ ${response} = "1" ]]; then
		met_file_type="processed_ll"
	    elif [[ ${response} = "2" ]]; then
		met_file_type="raw_ll"
		met_desc="raw"
	    else
		valid_response=0
		printf "Invalid option. Try again.\n"
	    fi
	done

	# Ask user to specify meteorology for advection.
	printf "${thinline}Choose meteorology for advection:${thinline}"
	printf "  1. 0.25x0.3125 3-hourly winds (Recommended)\n"
	printf "  2. C720 1-hourly derived winds (Recommended for stretched grid)\n"
	printf "  3. C720 1-hourly mass fluxes (Beta)\n"
	valid_response=0
	while [ "${valid_response}" -eq 0 ]; do
	    valid_response=1
	    read -p "${USER_PROMPT}" response
	    if [[ ${response} = "1" ]]; then
		adv_flux_src="3hr_wind"
	    elif [[ ${response} = "2" ]]; then
		adv_flux_src="1hr_derived_wind"
	    elif [[ ${response} = "3" ]]; then
		adv_flux_src="1hr_mass_flux"
	    else
		valid_response=0
		printf "Invalid option. Try again.\n"
	    fi
	done
	    
	# Set ExtData.rc settings for met data used for advection
	if [[ ${adv_flux_src} = "3hr_wind" ]]; then
	    if [[ ${met_file_type} = "processed_ll" ]]; then
		RUNDIR_VARS+="$(cat ${metSettingsDir}/geosfp/advection_met/geosfp.preprocessed_3hr_0.25x0.625_wind.txt)\n"
		RUNDIR_VARS+="$(cat ${metSettingsDir}/geosfp/advection_met/geosfp.preprocessed_3hr_0.25x0.625_PS_SPHU.txt)\n"
	    elif [[ ${met_file_type} = "raw_ll" ]]; then
		RUNDIR_VARS+="$(cat ${metSettingsDir}/geosfp/advection_met/geosfp.raw_3hr_0.25x0.625_wind.txt)\n"
		RUNDIR_VARS+="$(cat ${metSettingsDir}/geosfp/advection_met/geosfp.raw_3hr_0.25x0.625_PS_SPHU.txt)\n"
	    fi

	elif [[ ${adv_flux_src} = "1hr_derived_wind" ]]; then
	    RUNDIR_VARS+="$(cat ${metSettingsDir}/geosfp/advection_met/geosfp.derived_1hr_c720_wind_PS_SPHU.txt)\n"
	    
	elif [[ ${adv_flux_src} = "1hr_mass_flux" ]]; then
	    if [[ ${met_file_type} = "processed_ll" ]]; then
		RUNDIR_VARS+="$(cat ${metSettingsDir}/geosfp/advection_met/geosfp.raw_1hr_c720_mass_flux_PS_SPHU_C.txt)\n"
		RUNDIR_VARS+="$(cat ${metSettingsDir}/geosfp/advection_met/geosfp.preprocessed_3hr_0.25x0.625_wind.txt)\n"
	    elif [[ ${met_file_type} = "raw_ll" ]]; then
		RUNDIR_VARS+="$(cat ${metSettingsDir}/geosfp/advection_met/geosfp.raw_1hr_c720_mass_flux_PS_SPHU_C.txt)\n"
		RUNDIR_VARS+="$(cat ${metSettingsDir}/geosfp/advection_met/geosfp.raw_3hr_0.25x0.625_wind.txt)\n"
	    fi

	fi

	# Set ExtData.rc settings for met-fields used in GEOS-Chem
	if [[ ${met_file_type} = "raw_ll" ]]; then
	    RUNDIR_VARS+="$(cat ${metSettingsDir}/geosfp/geosfp.raw_ll.txt)\n"
	elif [[ ${met_file_type} = "processed_ll" ]]; then
	    RUNDIR_VARS+="$(cat ${metSettingsDir}/geosfp/geosfp.preprocessed_ll.txt)\n"
	fi
	
    elif [[ ${met_num} = "3" ]]; then
	
	met="geosit"
	
 	# Ask user to specify processed or raw files
	printf "${thinline}Choose meteorology files:${thinline}"
	printf "  1. C180 processed files from the GEOS-Chem data archive (Recommended)\n"
	printf "  2. C180 raw files downloaded from NASA GMAO\n"
	printf "  3. 0.5x0.625 processed files from the GEOS-Chem data archive\n"
	printf "  4. 0.5x0.625 raw files downloaded from NASA GMAO\n"
	printf "\nNOTES:\n"
	printf " - Processed files have unused variables removed, are vertically flipped,\n"
	printf "   and are concatenated into daily files. Using them speeds up GCHP runs.\n"
	printf " - Using input data on a cubed-sphere grid reduces regridding error in GCHP.\n\n"
	valid_response=0
	while [ "${valid_response}" -eq 0 ]; do
	    valid_response=1
	    read -p "${USER_PROMPT}" response
	    if [[ ${response} = "1" ]]; then
		met_file_type="processed_cs"
		met_desc="processed_cs"
	    elif [[ ${response} = "2" ]]; then
		met_file_type="raw_cs"
		met_desc="raw_cs"
	    elif [[ ${response} = "3" ]]; then
		met_file_type="processed_ll"
		met_desc="processed_ll"
	    elif [[ ${response} = "4" ]]; then
		met_file_type="raw_ll"
		met_desc="raw_ll"
	    else
		valid_response=0
		printf "Invalid option. Try again.\n"
	    fi
	done

	# If using cubed-sphere raw files ask user to specify meteoerology for
	# advection. If using raw lat-lon then always use winds.
	if [[ ${met_file_type} = "raw_ll" || ${met_file_type} = "processed_ll" ]]; then
	    adv_flux_src="3hr_wind"
	elif [[ ${met_file_type} = "processed_cs" || ${met_file_type} = "raw_cs" ]]; then
	    printf "${thinline}Choose meteorology for advection:${thinline}"
	    printf "  1. C180 3-hourly winds (Recommended)\n"
	    printf "  2. C180 1-hourly mass fluxes (Beta)\n"
	    valid_response=0
	    while [ "${valid_response}" -eq 0 ]; do
	        valid_response=1
	        read -p "${USER_PROMPT}" response
	        if [[ ${response} = "1" ]]; then
		    adv_flux_src="3hr_wind"
	        elif [[ ${response} = "2" ]]; then
		    adv_flux_src="1hr_mass_flux"
	        else
	    	    valid_response=0
	    	    printf "Invalid option. Try again.\n"
	        fi
	    done
	fi
	
	# If using raw files, ask user if they are using discover
	if [[ ${met_file_type} = "raw_cs" || ${met_file_type} = "raw_ll" ]]; then
	    printf "${thinline}Are you running on the NASA discover cluster? (y/n)${thinline}"
	    valid_response=0
	    while [ "$valid_response" -eq 0 ]; do
		read -p "${USER_PROMPT}" use_discover
		if [[ ${use_discover} = "y" ]]; then
		    valid_response=1
		elif [[ ${use_discover} = "n" ]]; then
		    valid_response=1
		else
		    printf "Invalid option. Try again.\n"
		fi
	    done
	else
	    use_discover=n
	fi
	
	# Set text files containing settings for met data. Different settings based options above.
	if [[ ${met_file_type} = "processed_ll" ]]; then
	    RUNDIR_VARS+="$(cat ${metSettingsDir}/geosit/geosit.preprocessed_0.5x0.625.txt)\n"
	
	elif [[ ${met_file_type} = "processed_cs" ]]; then

	    if [[ ${adv_flux_src} = "1hr_mass_flux" ]]; then
	        RUNDIR_VARS+="$(cat ${metSettingsDir}/geosit/advection_met/geosit.preprocessed_1hr_c180_mass_flux.txt)\n"
	    elif [[ ${adv_flux_src} == "3hr_wind" ]]; then
		RUNDIR_VARS+="$(cat ${metSettingsDir}/geosit/advection_met/geosit.preprocessed_3hr_c180_wind.txt)\n"
	    fi
	    RUNDIR_VARS+="$(cat ${metSettingsDir}/geosit/geosit.preprocessed_c180.txt)\n"

	else

	    # Using raw GEOS-IT data downloaded from GMAO
	    if [[ ${use_discover} = "y" ]]; then

		if [[ ${met_file_type} = "raw_ll" ]]; then
		    RUNDIR_VARS+="$(cat ${metSettingsDir}/geosit/discover/geosit.raw_0.5x0.625.txt)\n"
		else
	            if [[ ${adv_flux_src} = "1hr_mass_flux" ]]; then
			RUNDIR_VARS+="$(cat ${metSettingsDir}/geosit/discover/geosit.raw_1hr_c180_mass_flux.txt)\n"
	            elif [[ ${adv_flux_src} = "3hr_wind" ]]; then
			RUNDIR_VARS+="$(cat ${metSettingsDir}/geosit/discover/geosit.raw_3hr_c180_wind.txt)\n"
		    fi
		    RUNDIR_VARS+="$(cat ${metSettingsDir}/geosit/discover/geosit.raw_c180.txt)\n"
		fi

	    else
		# Not using discover
		if [[ ${met_file_type} = "raw_ll" ]]; then
		    RUNDIR_VARS+="$(cat ${metSettingsDir}/geosit/geosit.raw_0.5x0.625.txt)\n"
		else
	            if [[ ${adv_flux_src} = "1hr_mass_flux" ]]; then
			RUNDIR_VARS+="$(cat ${metSettingsDir}/geosit/advection_met/geosit.raw_1hr_c180_mass_flux.txt)\n"
	            elif [[ ${adv_flux_src} = "3hr_wind" ]]; then
			RUNDIR_VARS+="$(cat ${metSettingsDir}/geosit/advection_met/geosit.raw_3hr_c180_wind.txt)\n"
		    fi
		    RUNDIR_VARS+="$(cat ${metSettingsDir}/geosit/geosit.raw_c180.txt)\n"
		fi

	    fi
	fi  # end GEOS-IT
	
    else
	valid_met=0
	printf "Invalid meteorology option. Try again.\n"
    fi
done

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

    # Replace ~ with the user's home directory
    # NOTE: This is a safe algorithm.
    if [[ "${rundir_path}" =~ '~' ]]; then
	rundir_path="${rundir_path/#\~/$HOME}"
	echo "Expanding to: ${rundir_path}"
    fi

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
	rundir_name=gchp_${met}_${sim_name}
	if [[ "${sim_extra_option}" != "none" ]]; then
	    rundir_name=${rundir_name}_${sim_extra_option}
	fi
	if [[ "x${met_desc}" != "x" ]]; then
	    rundir_name=${rundir_name}_${met_desc}
	fi
	if [[ "x${adv_flux_src}" != "x" ]]; then
	    rundir_name=${rundir_name}_using_${adv_flux_src}
	fi
	printf "  -- Using default directory name ${rundir_name}\n"
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
cp ${gcdir}/run/shared/cleanRunDir.sh ${rundir}
cp ./archiveRun.sh                    ${rundir}
cp ./logging.yml                      ${rundir}
cp ./README.md                        ${rundir}
cp ./setEnvironmentLink.sh            ${rundir}
cp ./setRestartLink.sh                ${rundir}
cp ./checkRunSettings.sh              ${rundir}
cp ./gitignore                        ${rundir}/.gitignore

# Copy file to auto-update common settings. Use adjoint version for CO2.
cp ./setCommonRunSettings.sh.template  ${rundir}/setCommonRunSettings.sh

if [[ "x${sim_name}" == "xfullchem" || "x${sim_name}" == "xcarbon" ]]; then
    cp -r ${gcdir}/run/shared/metrics.py  ${rundir}
    chmod 744 ${rundir}/metrics.py
fi

# Copy the KPP standalone interface config file to ther rundir (fullchem only)
if [[ "x${sim_name}" == "xfullchem"  ]]; then
    cp -r ${gcdir}/run/shared/kpp_standalone_interface.yml ${rundir}
    chmod 644 ${rundir}/kpp_standalone_interface.yml
fi

# Set permissions
chmod 744 ${rundir}/cleanRunDir.sh
chmod 744 ${rundir}/archiveRun.sh
chmod 744 ${rundir}/setEnvironmentLink.sh
chmod 744 ${rundir}/setRestartLink.sh
chmod 744 ${rundir}/setCommonRunSettings.sh
chmod 744 ${rundir}/checkRunSettings.sh

# Copy species database; append APM or TOMAS species if needed
# Also copy APM input files to the run directory
cp -r ${gcdir}/run/shared/species_database.yml   ${rundir}
if [[ ${sim_extra_option} =~ "TOMAS" ]]; then
    cat ${gcdir}/run/shared/species_database_tomas.yml >> ${rundir}/species_database.yml
elif [[ ${sim_extra_option} =~ "APM" ]]; then
    cat ${gcdir}/run/shared/species_database_apm.yml >> ${rundir}/species_database.yml
    cp ${gcdir}/run/shared/apm_tmp.dat ${rundir}/apm_tmp.dat
    cp ${gcdir}/run/shared/input.apm   ${rundir}/input.apm
fi

# Create symbolic link to code directory
ln -s ${wrapperdir} ${rundir}/CodeDir
ln -s ${wrapperdir}/run/runScriptSamples ${rundir}/runScriptSamples

#--------------------------------------------------------------------
# Link to initial restart files, set start in cap_restart
#--------------------------------------------------------------------
restarts=${GC_DATA_ROOT}/GEOSCHEM_RESTARTS
if [[ "x${sim_name}" == "xfullchem" ]]; then
    if [[ ${sim_extra_option} =~ "TOMAS" ]]; then
	start_date='20190701'
	restart_dir='v2024-01'
	restart_name="${sim_extra_option}"
    else
	start_date='20190701'
	restart_dir='GC_14.5.0'
	restart_name="${sim_name}"
    fi
elif [[ "x${sim_name}" == "xtagO3" ]]; then
    # NOTE: we use the fullchem restart file for tagO3
    start_date='20190701'
    restart_dir='GC_14.5.0'
    restart_name="fullchem"
elif [[ "x${sim_name}" == "xTransportTracers" ]]; then
    start_date='20190101'
    restart_dir='GC_14.2.0'
    restart_name="${sim_name}"
elif [[ ${sim_name} = "carbon" ]]; then
    start_date='20190101'
    restart_dir='v2025-07'
    restart_name="${sim_name}"
fi
for N in 24 30 48 90 180
do
    # Do not include c24 and c48 if using GEOS-IT mass fluxes. MAPL cannot regrid
    # C180 mass fluxes to those grid resolutions.
    if [[ "${met}" == "geosit" && "${adv_flux_src}" == "mass_flux" ]]; then
	if [[ "$N" == "24" || "$N" == "48" ]]; then
	    continue
	fi
    fi
    old_prefix="GEOSChem.Restart.${restart_name}"
    new_prefix="GEOSChem.Restart"
    echo "${start_date} 000000" > ${rundir}/cap_restart
    initial_rst="${restarts}/${restart_dir}/${old_prefix}.${start_date}_0000z.c${N}.nc4"
    linkname="${rundir}/Restarts/${new_prefix}.${start_date}_0000z.c${N}.nc4"
    ln -s ${initial_rst} ${linkname}
done

#--------------------------------------------------------------------
# Navigate to run directory and set up input files
#--------------------------------------------------------------------
cd ${rundir}

RUNDIR_VARS+="RUNDIR_SIM_DUR_YYYYMMDD='00000100'\n"
RUNDIR_VARS+="RUNDIR_SIM_DUR_HHmmSS='000000'\n"

# Use monthly diagnostics by default
RUNDIR_VARS+="RUNDIR_HIST_TIME_AVG_DUR='010000'\n"
RUNDIR_VARS+="RUNDIR_HIST_TIME_AVG_FREQ='010000'\n"
RUNDIR_VARS+="RUNDIR_HIST_INST_DUR='010000'\n"
RUNDIR_VARS+="RUNDIR_HIST_INST_FREQ='010000'\n"
RUNDIR_VARS+="RUNDIR_HIST_MONTHLY_DIAG='1'\n"

# Set default compute resources
if [[ "${met}" == "geosit" && "${adv_flux_src}" == "1hr_mass_flux" ]]; then
    RUNDIR_VARS+="RUNDIR_NUM_CORES='24'\n"
    RUNDIR_VARS+="RUNDIR_NUM_NODES='1'\n"
    RUNDIR_VARS+="RUNDIR_CORES_PER_NODE='24'\n"
else
    RUNDIR_VARS+="RUNDIR_NUM_CORES='96'\n"
    RUNDIR_VARS+="RUNDIR_NUM_NODES='2'\n"
    RUNDIR_VARS+="RUNDIR_CORES_PER_NODE='48'\n"
fi

# Set default grid resolution
if [[ "${met}" == "geosit" && "${adv_flux_src}" == "1hr_mass_flux" ]]; then
    RUNDIR_VARS+="RUNDIR_CS_RES='30'\n"
else
    RUNDIR_VARS+="RUNDIR_CS_RES='24'\n"
fi

# Assign appropriate file paths and settings in HEMCO_Config.rc
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
RUNDIR_VARS+="$(cat ${metSettingsDir}/gmao_hemco.txt)\n"
if [[ "x${sim_extra_option}" == "xbenchmark"        ||
      "x${sim_extra_option}" == "xaciduptake"       ||
      "x${sim_extra_option}" == "xmarinePOA"        ||
      "x${sim_extra_option}" == "xcomplexSOA_SVPOA" ||
      "x${sim_extra_option}" == "xAPM"              ||
      "x${sim_name}"         == "xcarbon"           ||
      "x${sim_extra_option}" == "xTOMAS15"          ||
      "x${sim_extra_option}" == "xTOMAS40"          ||
      "x${sim_name}"         == "xPOPs"             ||
      "x${sim_name}"         == "xTransportTracers" ||
      "x${sim_name}"         == "xtagO3"        ]]; then
    RUNDIR_VARS+="RUNDIR_INITIAL_RESTART_SPECIES_REQUIRED='0'\n"
else
    RUNDIR_VARS+="RUNDIR_INITIAL_RESTART_SPECIES_REQUIRED='1'\n"
fi

#--------------------------------------------------------------------
# Replace settings in config files with RUNDIR variables
#--------------------------------------------------------------------

# Define a subdirectory for rundir configuration files
rundir_config_dirname=CreateRunDirLogs
rundir_config=${rundir}/${rundir_config_dirname}
mkdir -p ${rundir_config}

# Save RUNDIR variables to a file in the rundirConfig folder
rundir_config_logname=rundir_vars.txt
rundir_config_log=${rundir}/${rundir_config_dirname}/${rundir_config_logname}
echo -e "$RUNDIR_VARS" > ${rundir_config_log}

# Initialize run directory
${srcrundir}/init_rd.sh ${rundir_config_log}

# Call function to setup configuration files with settings common between
# GEOS-Chem Classic and GCHP. This script mainly now adds species to 
# geoschem_config.yml and modifies diagnostic output based on simulation type.
if [[ "x${sim_name}" = "xfullchem" ]]; then
    set_common_settings "${sim_extra_option}" "GCHP"
fi

# If necessary, edit config files for a carbon single species simulation
if [[ "x${sim_name}" == "xcarbon" ]]; then
    if [[ "x${sim_extra_option}" != "xnone" ]]; then
	singleCarbonSpecies "${sim_extra_option}" "${rundir}"
    fi
fi

# Call setCommonRunSettings.sh so that all config files are consistent with its
# default settings. Suppress informational prints.
chmod +x setCommonRunSettings.sh
./setCommonRunSettings.sh --silent

#--------------------------------------------------------------------
# Navigate back to source code directory
#--------------------------------------------------------------------
cd ${srcrundir}

#----------------------------------------------------------------------
# Archive repository version in run directory file rundir.version
#----------------------------------------------------------------------
version_log=${rundir}/${rundir_config_dirname}/rundir.version
echo "This run directory was created with ${srcrundir}/createRunDir.sh." > ${version_log}
echo " " >> ${version_log}
echo "GEOS-Chem repository version information:" >> ${version_log}
cd ${gcdir}
remote_url=$(git config --get remote.origin.url)
code_branch=$(git rev-parse --abbrev-ref HEAD)
last_commit=$(git log -n 1 --pretty=format:"%s")
commit_date=$(git log -n 1 --pretty=format:"%cd")
commit_user=$(git log -n 1 --pretty=format:"%cn")
commit_hash=$(git log -n 1 --pretty=format:"%h")
cd ${srcrundir}
printf "\n  Remote URL: ${remote_url}" >> ${version_log}
printf "\n  Branch: ${code_branch}"    >> ${version_log}
printf "\n  Commit: ${last_commit}"    >> ${version_log}
printf "\n  Date: ${commit_date}"      >> ${version_log}
printf "\n  User: ${commit_user}"      >> ${version_log}
printf "\n  Hash: ${commit_hash}"      >> ${version_log}

#-----------------------------------------------------------------
# Ask user whether to track run directory changes with git
#-----------------------------------------------------------------
printf "${thinline}Do you want to track run directory changes with git? (y/n)${thinline}"
valid_response=0
while [ "$valid_response" -eq 0 ]; do
    read -p "${USER_PROMPT}" enable_git
    if [[ ${enable_git} = "y" ]]; then
	cd ${rundir}
	printf "\n\nChanges to the following run directory files are tracked by git:\n\n" >> ${version_log}
	printf "\n"
	git init
	git add *.rc *.sh *.yml input.nml
	if [[ "x${sim_name}" == "xfullchem" || "x${sim_name}" == "xcarbon" ]]; then
	    git add *.py
	fi
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

printf "\n${thinline}Created ${rundir}\n"
printf "\n  -- See ${rundir_config_dirname}/${rundir_config_logname} for summary of default run directory settings"
printf "\n  -- This run directory is set up for simulation start date $start_date"
printf "\n  -- Restart files for this date at different grid resolutions are in the"
printf "\n     Restarts subdirectory"
printf "\n  -- To update start time, edit configuration file cap_restart and"
printf "\n     add or symlink file Restarts/GEOSChem.Restart.YYYYMMDD_HHmmz.cN.nc"
printf "\n     where YYYYMMDD_HHmm is start date and time"
printf "\n  -- Edit other commonly changed run settings in setCommonRunSettings.sh"
printf "\n  -- See build/README for compilation instructions"
printf "\n  -- Example run scripts are in the runScriptSamples subdirectory"
printf "\n  -- For more information visit the GCHP user guide at"
printf "\n     https://readthedocs.org/projects/gchp/\n\n"

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
	git commit -m "Update header of ${rundir_config_dirname}/${rundir_config_lognbame}" > /dev/null
	cd ${srcrundir}
    fi
fi

#-----------------------------------------------------------------
# Create and populate warnings file
#-----------------------------------------------------------------
if [[ $met == "geosfp" ]]; then
   echo -e ${fp_msg} > ${rundir}/warnings.txt
fi

exit 0
