#!/bin/bash

# createRunDir.sh: Create GCHP run directory
#
# Optional argument: run directory name
#
# If optional run directory name argument is not passed then the user
# will be prompted to enter a name interactively, or choose to use the
# default name gchp_{met}_{sim_name}.
#
# Usage: ./createRunDir.sh [rundirname]
#
# Initial version: E. Lundgren,10/5/2018

srcrundir=$(pwd -P)
cd ${srcrundir}
cd ../..
gcdir=$(pwd -P)
cd ../../../..
wrapperdir=$(pwd -P)
cd ${srcrundir}

# Load file with utility functions to setup configuration files
. ${gcdir}/run/shared/setupConfigFiles.sh

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
	read -e extdata
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

#-----------------------------------------------------------------
# Ask user to select simulation type
#-----------------------------------------------------------------
printf "${thinline}Choose simulation type:${thinline}"
printf "   1. Full chemistry\n"
printf "   2. TransportTracers\n"
printf "   3. CO2 w/ CMS-Flux emissions\n"
valid_sim=0
while [ "${valid_sim}" -eq 0 ]; do
    read sim_num
    valid_sim=1
    if [[ ${sim_num} = "1" ]]; then
	sim_name=fullchem
    elif [[ ${sim_num} = "2" ]]; then
	sim_name=TransportTracers
    elif [[ ${sim_num} = "3" ]]; then
	sim_name=CO2
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
	read sim_option
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
		read soa_option
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
		read tomas_option
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
            printf "*** IMPORTANT: You must manually specify -DRRTMG=y when compiling the model. ***\n"
	else
	    valid_sim_option=0
	    printf "Invalid simulation option. Try again.\n"
	fi
    done


# Currently no transport tracer extra options
elif [[ ${sim_name} = "TransportTracers" ]]; then
    sim_extra_option=none
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
    RUNDIR_VARS+="RUNDIR_VOLC_CLIMATOLOGY='\$ROOT/VOLCANO/v2021-09/so2_volcanic_emissions_CARN_v202005.degassing_only.rc'\n"

    if [[ "x${sim_extra_option}" == "xbenchmark" ]]; then
	RUNDIR_VARS+="RUNDIR_VOLC_TABLE='\$ROOT/VOLCANO/v2021-09/so2_volcanic_emissions_CARN_v202005.degassing_only.rc'\n"
    else
	RUNDIR_VARS+="RUNDIR_VOLC_TABLE='\$ROOT/VOLCANO/v2021-09/\$YYYY/\$MM/so2_volcanic_emissions_Carns.$YYYY$MM$DD.rc'\n"
    fi
fi

#-----------------------------------------------------------------
# Ask user to select meteorology source
#-----------------------------------------------------------------
printf "${thinline}Choose meteorology source:${thinline}"
printf "  1. MERRA-2 (Recommended)\n"
printf "  2. GEOS-FP \n"

valid_met=0
while [ "${valid_met}" -eq 0 ]; do
    read met_num
    valid_met=1
    if [[ ${met_num} = "1" ]]; then
	met="merra2"
	RUNDIR_VARS+="$(cat ${gcdir}/run/shared/settings/merra2.txt)\n"
    elif [[ ${met_num} = "2" ]]; then
	met="geosfp"
	RUNDIR_VARS+="$(cat ${gcdir}/run/shared/settings/geosfp.txt)\n"
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
    read -e rundir_path

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
            read mk_rundir
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
    read -e rundir_name
    if [[ -z "${rundir_name}" ]]; then
	if [[ "${sim_extra_option}" = "none" ]]; then
	    rundir_name=gchp_${met}_${sim_name}
	else
	    rundir_name=gchp_${met}_${sim_name}_${sim_extra_option}
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
	read -e new_rundir
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

# Define a subdirectory for rundir configuration files
rundir_config=${rundir}/CreateRunDirLogs
mkdir -p ${rundir_config}

# Copy run directory files and subdirectories
cp ${gcdir}/run/shared/cleanRunDir.sh ${rundir}
cp ./archiveRun.sh                    ${rundir}
cp ./logging.yml                      ${rundir}
cp ./README                           ${rundir}
cp ./setEnvironmentLink.sh            ${rundir}
cp ./setRestartLink.sh                ${rundir}
cp ./checkRunSettings.sh              ${rundir}
cp ./gitignore                        ${rundir}/.gitignore

# Copy file to auto-update common settings. Use adjoint version for CO2.
cp ./setCommonRunSettings.sh.template  ${rundir}/setCommonRunSettings.sh

if [[ "x${sim_name}" == "xfullchem" || "x${sim_name}" == "xCH4" ]]; then
    cp -r ${gcdir}/run/shared/metrics.py  ${rundir}
    chmod 744 ${rundir}/metrics.py
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
for N in 24 48 90 180 360
do
    old_prefix="GCHP.Restart.${sim_name}"
    new_prefix="GEOSChem.Restart"
    if [[ ${sim_name} = "fullchem" ]]; then
        start_date='20190701'
        restart_dir='v2021-09'
    elif [[ ${sim_name} = "TransportTracers" ]]; then
        start_date='20190101'
        restart_dir='GC13.0.0'
    fi
    echo "${start_date} 000000" > ${rundir}/cap_restart
    initial_rst="${restarts}/${restart_dir}/${old_prefix}.${start_date}_0000z.c${N}.nc4"
    linkname="${rundir}/Restarts/${new_prefix}.${start_date}_000000z.c${N}.nc4"
    ln -s ${initial_rst} ${linkname}
done

# Add restart file to RUNDIR vars (need to fix/test this)
RUNDIR_VARS+="RUNDIR_RESTART_FILE=${new_prefix}.${start_date}.c"'{CS_RES}'".nc4\n"

#--------------------------------------------------------------------
# Navigate to run directory and set up input files
#--------------------------------------------------------------------
cd ${rundir}

# Special handling for start/end date based on simulation so that
# start year/month/day matches default initial restart file.
if [[ "x${sim_name}" == "xTransportTracers" ]]; then
    startdate='20190101'
    enddate='20190201'
elif [[ "x${sim_name}" == "xCO2" ]]; then
    startdate='20140901'
    enddate='20141001'
else
    startdate='20190701'
    enddate='20190801'
fi
RUNDIR_VARS+="RUNDIR_SIM_START_DATE=$startdate\n"
RUNDIR_VARS+="RUNDIR_SIM_END_DATE=$enddate\n"
RUNDIR_VARS+="RUNDIR_SIM_START_TIME='000000'\n"
RUNDIR_VARS+="RUNDIR_SIM_END_TIME='000000'\n"
RUNDIR_VARS+="RUNDIR_SIM_DUR_YYYYMMDD='00000100'\n"
RUNDIR_VARS+="RUNDIR_SIM_DUR_HHmmSS='000000'\n"

# Use monthly diagnostics by default
RUNDIR_VARS+="RUNDIR_HIST_TIME_AVG_DUR='7440000'\n"
RUNDIR_VARS+="RUNDIR_HIST_TIME_AVG_FREQ='7440000'\n"
RUNDIR_VARS+="RUNDIR_HIST_INST_DUR='7440000'\n"
RUNDIR_VARS+="RUNDIR_HIST_INST_FREQ='7440000'\n"
RUNDIR_VARS+="RUNDIR_HIST_MONTHLY_DIAG='1'\n"

# Set default compute resources
RUNDIR_VARS+="RUNDIR_NUM_CORES='96'\n"
RUNDIR_VARS+="RUNDIR_NUM_NODES='2'\n"
RUNDIR_VARS+="RUNDIR_CORES_PER_NODE='48'\n"
RUNDIR_VARS+="RUNDIR_CS_RES='24'\n"

# Turn on GEOS-Chem timers for benchmark simulations
if [[ "${sim_extra_option}" == "benchmark" ]]; then
    RUNDIR_VARS+="RUNDIR_USE_GCCLASSIC_TIMERS='true'\n"
else
    RUNDIR_VARS+="RUNDIR_USE_GCCLASSIC_TIMERS='false'\n"
fi

# Assign appropriate file paths and settings in HEMCO_Config.rc
if [[ "${sim_extra_option}" == "benchmark" ]]; then
    RUNDIR_VARS+="RUNDIR_DUSTDEAD_EXT='on '\n"
    RUNDIR_VARS+="RUNDIR_SEASALT_EXT='on '\n"
    RUNDIR_VARS+="RUNDIR_SOILNOX_EXT='on '\n"
    RUNDIR_VARS+="RUNDIR_OFFLINE_DUST='false'\n"
    RUNDIR_VARS+="RUNDIR_OFFLINE_BIOVOC='false'\n"
    RUNDIR_VARS+="RUNDIR_OFFLINE_SEASALT='false'\n"
    RUNDIR_VARS+="RUNDIR_OFFLINE_SOILNOX='false'\n"
else
    if [[ "${sim_extra_option}" == "marinePOA" ]]; then
	RUNDIR_VARS+="RUNDIR_SEASALT_EXT='on '\n"
	RUNDIR_VARS+="RUNDIR_OFFLINE_SEASALT='false'\n"
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
    RUNDIR_VARS+="RUNDIR_SOILNOX_EXT='off'\n"
    RUNDIR_VARS+="RUNDIR_OFFLINE_BIOVOC='true '\n"
    RUNDIR_VARS+="RUNDIR_OFFLINE_SOILNOX='true '\n"
fi
RUNDIR_VARS+="$(cat ${gcdir}/run/shared/settings/gmao_hemco.txt)\n"

#--------------------------------------------------------------------
# Replace settings in config files with RUNDIR variables
#--------------------------------------------------------------------

# Save RUNDIR variables to file
rundir_config_log=${rundir_config}/rundir_vars.txt
echo -e "$RUNDIR_VARS" > ${rundir_config_log}
sort -o ${rundir_config_log} ${rundir_config_log}

# Initialize run directory
${srcrundir}/init_rd.sh ${rundir_config_log}

#--------------------------------------------------------------------
# Print run direcory setup info to screen
#--------------------------------------------------------------------

printf "\n  -- This run directory has been set up for $startdate - $enddate.\n"
printf "\n  -- The default frequency and duration of diagnostics is set to monthly.\n"
printf "\n  -- You may modify these settings in setCommonRunSettings.sh.\n"

# Call function to setup configuration files with settings common between
# GEOS-Chem Classic and GCHP. This script mainly now adds species to 
# geoschem_config.yml and modifies diagnostic output based on simulation type.
if [[ "x${sim_name}" = "xfullchem" ]]; then
    set_common_settings ${sim_extra_option}
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
version_log=${rundir_config}/rundir.version
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
    read enable_git
    if [[ ${enable_git} = "y" ]]; then
	cd ${rundir}
	printf "\n\nChanges to the following run directory files are tracked by git:\n\n" >> ${version_log}
	printf "\n"
	git init
	git add *.rc *.sh *.yml input.nml
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
# Done!
#-----------------------------------------------------------------
printf "\nCreated ${rundir}\n"

exit 0
