#!/bin/bash

# createRunDir.sh: Create GCHP run directory
#
# Optional argument: run directory name
#
# If optional run directory name argument is not passed then the user
# will be prompted to enter a name interactively, or choose to use the
# default name gchp_{simulation}.
#
# Usage: ./createRunDir.sh [rundirname]
#
# Initial version: E. Lundgren,10/5/2018

srcrundir=$(pwd -P)
cd ${srcrundir}
cd ../..
gcdir=$(pwd)
cd ../../../..
gchpdir=$(pwd)
cd ${srcrundir}

# Load file with utility functions to setup configuration files
. ${gcdir}/run/shared/setupConfigFiles.sh

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
    while [ "$valid_path" -eq 0 ]
    do
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

#-----------------------------------------------------------------
# Ask user to select simulation type
#-----------------------------------------------------------------
printf "${thinline}Choose simulation type:${thinline}"
printf "   1. Full chemistry\n"
printf "   2. TransportTracers\n"

valid_sim=0
while [ "${valid_sim}" -eq 0 ]; do
    read sim_num
    valid_sim=1
    if [[ ${sim_num} = "1" ]]; then
	sim_name=fullchem
    elif [[ ${sim_num} = "2" ]]; then
	sim_name=TransportTracers
    else
        valid_sim=0
	printf "Invalid simulation option. Try again.\n"
    fi
done

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
    printf "  8. Standard w/ RRTMG\n"
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
	    printf "  1. TOASS with 40 bins\n"
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
	met_name='MERRA2'
	met_name_lc="merra2"
	met_dir='MERRA2'
	met_resolution='05x0625'
	met_native='0.5x0.625'
	met_latres='05'
	met_lonres='0625'
	met_extension='nc4'
	met_cn_year='2015'
	pressure_unit='Pa '
	pressure_scale='0.01'
	dust_sf='3.86e-4'
    elif [[ ${met_num} = "2" ]]; then
	met_name='GEOSFP'
	met_name_lc="merra2"
	met_dir='GEOS_FP'
	met_resolution='025x03125'
	met_native='0.25x0.3125'
	met_latres='025'
	met_lonres='03125'
	met_extension='nc'
	met_cn_year='2011'
	pressure_unit='hPa'
	pressure_scale='1.0 '
	dust_sf='6.42e-5'
    else
	valid_met=0
	printf "Invalid meteorology option. Try again.\n"
    fi
done

#-----------------------------------------------------------------
# Ask user to define path where directoy will be created
#-----------------------------------------------------------------
printf "${thinline}Enter path where the run directory will be created:${thinline}"
valid_path=0
while [ "$valid_path" -eq 0 ]; do
    read -e rundir_path

    # Test for quitting
    if [[ "x${rundir_path}" = "xq" ]]; then
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
	    rundir_name=gchp_${met_name_lc}_${sim_name}
	else
	    rundir_name=gchp_${met_name_lc}_${sim_name}_${sim_extra_option}
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

# Copy run directory files and subdirectories
cp ${gcdir}/run/shared/cleanRunDir.sh ${rundir}
cp ./archiveRun.sh                    ${rundir}
cp ./input.nml                        ${rundir}
cp ./README                           ${rundir}
cp ./setEnvironment.sh                ${rundir}
cp ./gitignore                        ${rundir}/.gitignore
cp ./GCHP.rc.template                 ${rundir}/GCHP.rc
cp ./CAP.rc.template                  ${rundir}/CAP.rc
cp ./runConfig.sh.template            ${rundir}/runConfig.sh
cp ./input.geos.templates/input.geos.${sim_name}            ${rundir}/input.geos
cp ./HISTORY.rc.templates/HISTORY.rc.${sim_name}            ${rundir}/HISTORY.rc
cp ./ExtData.rc.templates/ExtData.rc.${sim_name}            ${rundir}/ExtData.rc
cp ./HEMCO_Config.rc.templates/HEMCO_Config.rc.${sim_name}  ${rundir}/HEMCO_Config.rc
cp ./HEMCO_Diagn.rc.templates/HEMCO_Diagn.rc.${sim_name}    ${rundir}/HEMCO_Diagn.rc
cp -r ./utils ${rundir}
if [[ ${sim_name} = "fullchem" ]]; then
    cp -r ${gcdir}/run/shared/metrics.py  ${rundir}
    chmod 744 ${rundir}/metrics.py
fi
mkdir ${rundir}/OutputDir

# Set permissions
chmod 744 ${rundir}/setEnvironment.sh
chmod 744 ${rundir}/cleanRunDir.sh
chmod 744 ${rundir}/runConfig.sh
chmod 744 ${rundir}/archiveRun.sh

# Copy species database; append APM or TOMAS species if needed
cp -r ${gcdir}/run/shared/species_database.yml   ${rundir}
if [[ ${sim_extra_option} =~ "TOMAS" ]]; then
    cat ${gcdir}/run/shared/species_database_tomas.yml >> ${rundir}/species_database.yml
elif [[ ${sim_extra_option} =~ "APM" ]]; then
    cat ${gcdir}/run/shared/species_database_apm.yml >> ${rundir}/species_database.yml
fi

# If benchmark simulation, put run script in directory
if [[ ${sim_extra_option} = "benchmark" ]]; then
    cp ${gcdir}/run/GCHP/runScriptSamples/operational_examples/harvard_gcst/gchp.benchmark.run ${rundir}
    chmod 744 ${rundir}/gchp.benchmark.run
fi

# Create symbolic links to data directories, restart files, code, run scripts
ln -s ${gchpdir}                                ${rundir}/CodeDir
ln -s ${gcdir}/run/GCHP/runScriptSamples        ${rundir}/runScriptSamples
ln -s ${GC_DATA_ROOT}/CHEM_INPUTS               ${rundir}/ChemDir
ln -s ${GC_DATA_ROOT}/HEMCO                     ${rundir}/HcoDir
if [ "${met_name}" == "GEOSFP" ]; then
   ln -s ${GC_DATA_ROOT}/GEOS_0.25x0.3125/GEOS_FP  ${rundir}/MetDir
else
   ln -s ${GC_DATA_ROOT}/GEOS_0.5x0.625/MERRA2  ${rundir}/MetDir
fi
restarts=${GC_DATA_ROOT}/GEOSCHEM_RESTARTS
for N in 24 48 90 180 360
do
    src_prefix="GCHP.Restart.${sim_name}."
    src_suffix=".c${N}.nc4"
    target_name=initial_GEOSChem_rst.c${N}_${sim_name}.nc
    if [[ ${sim_name} = "fullchem" ]]; then
        start_date="20190701_0000z"
        src_name="${src_prefix}${start_date}${src_suffix}"
        ln -s ${restarts}/GC_13.0.0/${src_name} ${rundir}/${target_name}
    elif [[ ${sim_name} = "TransportTracers" ]]; then
        start_date="20190101_0000z"
        src_name="${src_prefix}${start_date}${src_suffix}"
        ln -s ${restarts}/GC_13.0.0/${src_name} ${rundir}/${target_name}
    fi
done

#--------------------------------------------------------------------
# Navigate to run directory and set up input files
#--------------------------------------------------------------------
cd ${rundir}

# Replace token strings in certain files
sed -i -e "s|{SIMULATION}|${sim_name}|"       GCHP.rc
sed -i -e "s|{SIMULATION}|${sim_name}|"       runConfig.sh
sed -i -e "s|{DATA_ROOT}|${GC_DATA_ROOT}|"    input.geos
sed -i -e "s|{MET}|${met_name}|"              input.geos
sed -i -e "s|{SIM}|${sim_name}|"              input.geos
sed -i -e "s|{DATA_ROOT}|${GC_DATA_ROOT}|"    HEMCO_Config.rc
sed -i -e "s|{NATIVE_RES}|${met_native}|"     HEMCO_Config.rc
sed -i -e "s|{LATRES}|${met_latres}|"         HEMCO_Config.rc
sed -i -e "s|{LONRES}|${met_lonres}|"         HEMCO_Config.rc
sed -i -e "s|{DUST_SF}|${dust_sf}|"           HEMCO_Config.rc
sed -i -e "s|{MET_SOURCE}|${met_name}|"       ExtData.rc # 1st in line
sed -i -e "s|{MET_SOURCE}|${met_name}|"       ExtData.rc # 2nd in line
sed -i -e "s|{MET_RES}|${met_resolution}|"    ExtData.rc
sed -i -e "s|{NATIVE_RES}|${met_native}|"     ExtData.rc
sed -i -e "s|{LATRES}|${met_latres}|"         ExtData.rc
sed -i -e "s|{LONRES}|${met_lonres}|"         ExtData.rc
sed -i -e "s|{MET_EXT}|${met_extension}|"     ExtData.rc
sed -i -e "s|{MET_CN_YR}|${met_cn_year}|"     ExtData.rc # 1st in line
sed -i -e "s|{MET_CN_YR}|${met_cn_year}|"     ExtData.rc # 2nd in line
sed -i -e "s|{PRES_UNIT}|${pressure_unit}|"   ExtData.rc
sed -i -e "s|{PRES_SCALE}|${pressure_scale}|" ExtData.rc

# Special handling for start/end date based on simulation so that
# start year/month/day matches default initial restart file.
if [[ ${sim_extra_option} = "benchmark" ]]; then
    startdate="20190701"
    enddate="20190801"
elif [[ ${sim_name} = "fullchem" ]]; then
    startdate="20190701"
    enddate="20190701"
else
    startdate="20190101"
    enddate="20190201"
fi
sed -i -e "s|{DATE1}|${startdate}|"     ${rundir}/runConfig.sh
sed -i -e "s|{DATE2}|${enddate}|"       ${rundir}/runConfig.sh
sed -i -e "s|{DATE1}|${startdate}|"     ${rundir}/CAP.rc
sed -i -e "s|{DATE2}|${enddate}|"       ${rundir}/CAP.rc

# Special handling for benchmark simulation
if [[ ${sim_extra_option} = "benchmark" || ${sim_name} == "TransportTracers" ]]; then
    total_cores=48
    num_nodes=2
    num_cores_per_node=24
    grid_res=48
    timeAvg_freq="7440000"
    inst_freq="7440000"
    start_time="000000"
    end_time="000000"
    dYYYYMMDD="00000100"
    dHHmmSS="000000"
    printf "\n  -- This run directory has been set up for $startdate $start_time - $enddate $end_time."
    printf "\n  -- The default diagnostic frequency and duration is 31 days."
else
    total_cores=24
    num_nodes=1
    num_cores_per_node=24
    grid_res=24
    timeAvg_freq="010000"
    inst_freq="010000"
    start_time="000000"
    end_time="010000"
    dYYYYMMDD="00000000"
    dHHmmSS="010000"
    printf "\n  -- This run directory has been set up for $startdate $start_time - $enddate $end_time."
    printf "\n  -- The default diagnostic frequency and duration is hourly."
fi
printf "\n  -- You may modify these settings in runConfig.sh.\n"
timeAvg_dur=${timeAvg_freq}
inst_dur=${inst_freq}
sed -i -e "s|{TotalCores}|${total_cores}|"             ${rundir}/runConfig.sh
sed -i -e "s|{NumNodes}|${num_nodes}|"                 ${rundir}/runConfig.sh
sed -i -e "s|{NumCoresPerNode}|${num_cores_per_node}|" ${rundir}/runConfig.sh
sed -i -e "s|{GridRes}|${grid_res}|"                   ${rundir}/runConfig.sh
sed -i -e "s|{InstFreq}|${inst_freq}|"                 ${rundir}/runConfig.sh
sed -i -e "s|{InstDur}|${inst_dur}|"                   ${rundir}/runConfig.sh
sed -i -e "s|{AvgFreq}|${timeAvg_freq}|"               ${rundir}/runConfig.sh
sed -i -e "s|{AvgDur}|${timeAvg_dur}|"                 ${rundir}/runConfig.sh
sed -i -e "s|{TIME1}|${start_time}|"                   ${rundir}/runConfig.sh
sed -i -e "s|{TIME2}|${end_time}|"                     ${rundir}/runConfig.sh
sed -i -e "s|{dYYYYMMDD}|${dYYYYMMDD}|"                ${rundir}/runConfig.sh
sed -i -e "s|{dHHmmss}|${dHHmmSS}|"                    ${rundir}/runConfig.sh
sed -i -e "s|{TIME1}|${start_time}|"                   ${rundir}/CAP.rc
sed -i -e "s|{TIME2}|${end_time}|"                     ${rundir}/CAP.rc
sed -i -e "s|{dYYYYMMDD}|${dYYYYMMDD}|"                ${rundir}/CAP.rc
sed -i -e "s|{dHHmmss}|${dHHmmSS}|"                    ${rundir}/CAP.rc

# Call function to setup configuration files with settings common between
# GEOS-Chem Classic and GCHP.
if [[ ${sim_name} = "fullchem" ]]; then
    set_common_settings ${sim_extra_option}
fi

# Call runConfig.sh so that all config files are consistent with its
# default settings. Suppress informational prints.
./runConfig.sh --silent

#--------------------------------------------------------------------
# Navigate back to source code directory
#--------------------------------------------------------------------
cd ${srcrundir}

#----------------------------------------------------------------------
# Archive repository version in run directory file rundir.version
#----------------------------------------------------------------------
version_log=${rundir}/rundir.version
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
	git add *.rc *.sh *.yml input.geos input.nml
	if [[ ${sim_name} = "fullchem" ]]; then
            git add *.py
	fi
	git add README .gitignore
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
