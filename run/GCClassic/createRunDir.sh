#!/bin/bash

# createRunDir.sh: Create GEOS-Chem Classic run directory
#
# Optional argument: run directory name
#
# If optional run directory name argument is not passed then the user
# will be prompted to enter a name interactively, or choose to use the
# default name gchp_{simulation}/
#
# Usage: ./createRunDir.sh [rundirname]
#
# Initial version: M. Sulprizio, 6/24/2020

curdir=$(pwd)
cd ../..
gcdir=$(pwd)
cd ${curdir}

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
    printf "\nDefine path to ExtData."
    printf "\nThis will be stored in ${HOME}/.geoschem/config for future automatic use.\n"
    mkdir -p ${HOME}/.geoschem
fi

#-----------------------------------------------------------------
# One-time configuration of data root path in ~/.geoschem/config
#-----------------------------------------------------------------
if [[ -z "${GC_DATA_ROOT}" ]]; then
    printf "\nEnter path for ExtData:\n"
    valid_path=0
    while [ "$valid_path" -eq 0 ]
    do
	read extdata
	if [[ ${extdata} = "q" ]]; then
	    printf "\nExiting.\n"
	    exit 1
	elif [[ ! -d ${extdata} ]]; then
            printf "\nError: ${extdata} does not exist. Enter a new path or hit q to quit.\n"
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
printf "\nChoose simulation type:\n"
printf "   1. Standard\n"
printf "   2. Tropchem\n"
printf "   3. Complex SOA\n"
printf "   4. Complex SOA + SVPOA\n"
printf "   5. TOMAS15\n"
printf "   6. TOMAS40\n"
printf "   7. APM\n"
printf "   8. RRTMG\n"
printf "   9. Acid uptake\n"
printf "  10. Marine POA\n"
printf "  11. Benchmark\n"
printf "  12. Aerosol only\n"
printf "  13. CH4\n"
printf "  14. CO2\n"
printf "  15. Hg\n"
printf "  16. POPs: BaP\n"
printf "  17. POPs: PHE\n"
printf "  18. POPS: PYR\n"
printf "  19. tagCH4\n"
printf "  20. tagCO\n"
printf "  21. tagO3\n"
printf "  22. TransportTracers\n"

valid_sim=0
while [ "${valid_sim}" -eq 0 ]
do
    read sim_num
    if [[ ${sim_num} = "1" ]]; then
	sim_name=standard
	sim_type=${sim_name}
	valid_sim=1
    elif [[ ${sim_num} = "2" ]]; then
	sim_name=tropchem
	sim_type=${sim_name}
	valid_sim=1
    elif [[ ${sim_num} = "3" ]]; then
	sim_name=complexSOA
	sim_type=tropchem
	valid_sim=1
    elif [[ ${sim_num} = "4" ]]; then
	sim_name=complexSOA_SVPOA
	sim_type=tropchem
	valid_sim=1
    elif [[ ${sim_num} = "5" ]]; then
	sim_name=TOMAS15
	sim_type=tropchem
	valid_sim=1
    elif [[ ${sim_num} = "6" ]]; then
	sim_name=TOMAS40
	sim_type=tropchem
	valid_sim=1
    elif [[ ${sim_num} = "7" ]]; then
	sim_name=APM
	sim_type=tropchem
	valid_sim=1
    elif [[ ${sim_num} = "8" ]]; then
	sim_name=RRTMG
	sim_type=tropchem
	valid_sim=1
    elif [[ ${sim_num} = "9" ]]; then
	sim_name=aciduptake
	sim_type=tropchem
	valid_sim=1
    elif [[ ${sim_num} = "10" ]]; then
	sim_name=marinePOA
	sim_type=tropchem
	valid_sim=1
    elif [[ ${sim_num} = "11" ]]; then
	sim_name=benchmark
	sim_type=${sim_name}
	valid_sim=1
    elif [[ ${sim_num} = "12" ]]; then
	sim_name=aerosol
	sim_type=${sim_name}
	valid_sim=1
    elif [[ ${sim_num} = "13" ]]; then
	sim_name=CH4
	sim_type=${sim_name}
	valid_sim=1
    elif [[ ${sim_num} = "14" ]]; then
	sim_name=CO2
	sim_type=${sim_name}
	valid_sim=1
    elif [[ ${sim_num} = "15" ]]; then
	sim_name=Hg
	sim_type=${sim_name}
	valid_sim=1
    elif [[ ${sim_num} = "16" ]]; then
	sim_name=POPs_BaP
	sim_type=POPs
	valid_sim=1
    elif [[ ${sim_num} = "17" ]]; then
	sim_name=POPs_PHE
	sim_type=POPs
	valid_sim=1
    elif [[ ${sim_num} = "18" ]]; then
	sim_name=POPs_PYR
	sim_type=POPs
	valid_sim=1
    elif [[ ${sim_num} = "19" ]]; then
	sim_name=tagCH4
	sim_type=${sim_name}
	valid_sim=1
    elif [[ ${sim_num} = "20" ]]; then
	sim_name=tagCO
	sim_type=${sim_name}
	valid_sim=1
    elif [[ ${sim_num} = "21" ]]; then
	sim_name=tagO3
	sim_type=${sim_name}
	valid_sim=1
    elif [[ ${sim_num} = "22" ]]; then
	sim_name=TransportTracers
	sim_type=${sim_name}
	valid_sim=1
    else
	printf "Invalid simulation option. Try again.\n"
    fi
done

#-----------------------------------------------------------------
# Ask user to select meteorology source
#-----------------------------------------------------------------
printf "\nChoose meteorology source:\n"
printf "  1. MERRA2 (Recommended)\n"
printf "  2. GEOS-FP \n"
valid_met=0
while [ "${valid_met}" -eq 0 ]
do
    read met_num
    if [[ ${met_num} = "1" ]]; then
	met_name='MERRA2'
	met_dir='MERRA2'
	met_resolution='05x0625'
	met_native='0.5x0.625'
	met_latres='05'
	met_lonres='0625'
	met_extension='nc4'
	met_cn_year='2015'
	pressure_unit='Pa '
	pressure_scale='0.01'
	valid_met=1
	dust_sf='3.86e-4'
    elif [[ ${met_num} = "2" ]]; then
	met_name='GEOSFP'
	met_dir='GEOS_FP'
	met_resolution='025x03125'
	met_native='0.25x0.3125'
	met_latres='025'
	met_lonres='03125'
	met_extension='nc'
	met_cn_year='2011'
	pressure_unit='hPa'
	pressure_scale='1.0 '
	valid_met=1
	dust_sf='6.42e-5'
    else
	printf "Invalid meteorology option. Try again.\n"
    fi
done

#-----------------------------------------------------------------
# Ask user to select horizontal resolution
#-----------------------------------------------------------------
printf "\nChoose horizontal resolution:\n"
printf "  1. 4.0  x 5.0\n"
printf "  2. 2.0  x 2.5\n"
printf "  3. 0.5  x 0.625\n"
if [[ ${met_name} = "GEOSFP" ]]; then
    printf "  4. 0.25 x 0.3125\n"
fi

valid_res=0
while [ "${valid_res}" -eq 0 ]
do
    read res_num
    if [[ ${res_num} = "1" ]]; then
	grid_res='4x5'
        grid_res_long='4.0x5.0'
	valid_res=1
    elif [[ ${res_num} = "2" ]]; then
	grid_res='2x25'
	grid_res_long='2.0x2.5'
	valid_res=1
    elif [[ ${res_num} = "3" ]]; then
	grid_res='05x0625'
	grid_res_long='0.5x0.625'
	valid_res=1
    elif [[ ${res_num} = "4" ]]; then
	grid_res='025x03125'
	grid_res_long='0.25x0.3125'
	valid_res=1
    else
	printf "Invalid horizontal resolution option. Try again.\n"
    fi
done

#-----------------------------------------------------------------
# Ask user to select vertical resolution
#-----------------------------------------------------------------
printf "\nChoose number of levels:\n"
printf "  1. 72 (native)\n"
printf "  2. 47 (reduced)\n"

valid_lev=0
while [ "${valid_lev}" -eq 0 ]
do
    read lev_num
    if [[ ${lev_num} = "1" ]]; then
	grid_lev='72'
	valid_lev=1
    elif [[ ${res_num} = "2" ]]; then
	grid_lev='47'
	valid_lev=1
    else
	printf "Invalid vertical resolution option. Try again.\n"
    fi
done

#-----------------------------------------------------------------
# Ask user to define path where directoy will be created
#-----------------------------------------------------------------
printf "\nEnter path where the run directory will be created:\n"
valid_path=0
while [ "$valid_path" -eq 0 ]
do
    read rundir_path
    if [[ ${rundir_path} = "q" ]]; then
	printf "\nExiting.\n"
	exit 1
    elif [[ ! -d ${rundir_path} ]]; then
        printf "\nError: ${rundir_path} does not exist. Enter a new path or hit q to quit.\n"
    else
	valid_path=1
    fi
done

#-----------------------------------------------------------------
# Ask user to define run directoy name if not passed as argument
#-----------------------------------------------------------------
if [ -z "$1" ]; then
    printf "\nEnter run directory name, or press return to use default:\n"
    read rundir_name
    if [[ -z "${rundir_name}" ]]; then
	rundir_name=GC_${grid_res}_${met_name}_${sim_name}
	printf "Using default directory name ${rundir_name}\n"
    fi
else
    rundir_name=$1
fi

#-----------------------------------------------------------------
# Ask user for a new run directory name if specified one exists
#-----------------------------------------------------------------
rundir=${rundir_path}/${rundir_name}
valid_rundir=0
while [ "${valid_rundir}" -eq 0 ]
do
    if [[ -d ${rundir} ]]; then
	printf "Warning! ${rundir} already exists.\n"
        printf "Enter a different run directory name, or q to quit:\n"
	read new_rundir
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

#-----------------------------------------------------------------
# Copy run directory files and subdirectories
#-----------------------------------------------------------------
cp -r ./getRunInfo               ${rundir}
cp -r ../shared/download_data.py ${rundir}
cp -r ./runScriptSamples         ${rundir}
cp ./archiveRun.sh               ${rundir}
cp ./cleanRundir.sh              ${rundir}
cp ./setCodeDir.sh               ${rundir}
cp ./README                      ${rundir}
cp ./gitignore                   ${rundir}/.gitignore
cp ./input.geos.templates/input.geos.${sim_name}            ${rundir}/input.geos
cp ./HISTORY.rc.templates/HISTORY.rc.${sim_type}            ${rundir}/HISTORY.rc
cp ./HEMCO_Config.rc.templates/HEMCO_Config.rc.${sim_name}  ${rundir}/HEMCO_Config.rc
if [[ ${sim_type} =~ "chem" ]]; then
    cp ./HEMCO_Diagn.rc.templates/HEMCO_Diagn.rc.standard   ${rundir}/HEMCO_Diagn.rc
else
    cp ./HEMCO_Diagn.rc.templates/HEMCO_Diagn.rc.${sim_name}    ${rundir}/HEMCO_Diagn.rc
fi
mkdir ${rundir}/OutputDir

# Copy species database; append APM or TOMAS species if needed
cp -r ../shared/species_database.yml   ${rundir}
if [[ ${sim_name} =~ "TOMAS" ]]; then
    cat ../shared/species_database_tomas.yml >> ${rundir}/species_database.yml
elif [[ ${sim_name} =~ "APM" ]]; then
    cat ../shared/species_database_apm.yml >> ${rundir}/species_database.yml
fi

# If benchmark simulation, put run script in directory
if [ "${sim_name}" == "benchmark" ]; then
    cp ./runScriptSamples/geoschem.benchmark.run ${rundir}
    chmod 744 ${rundir}/geoschem.benchmark.run
fi

#--------------------------------------------------------------------
# Copy sample restart file
#--------------------------------------------------------------------
restarts=${GC_DATA_ROOT}/GEOSCHEM_RESTARTS/v2018-11
cp ${restarts}/initial_GEOSChem_rst.${grid_res}_${sim_name}.nc ${rundir}

#--------------------------------------------------------------------
# Create symbolic link to code directory
#--------------------------------------------------------------------
ln -s ${gcdir}                                  ${rundir}/CodeDir

#-----------------------------------------------------------------
# Replace token strings in certain files
#-----------------------------------------------------------------
sed -i -e "s|{DATA_ROOT}|${GC_DATA_ROOT}|"   ${rundir}/input.geos
sed -i -e "s|{MET}|${met_name}|"             ${rundir}/input.geos
sed -i -e "s|{SIM}|${sim_name}|"             ${rundir}/input.geos
sed -i -e "s|{RES}|${grid_res_long}|"        ${rundir}/input.geos
sed -i -e "s|{NLEV}|${grid_lev}|"            ${rundir}/input.geos
## TO DO: Add if nested... else here; for now default to global
sed -i -e "s|{LON_RANGE}|-180.0 180.0|"      ${rundir}/input.geos
sed -i -e "s|{LAT_RANGE}| -90.0  90.0|"      ${rundir}/input.geos
sed -i -e "s|{HALF_POLAR}|T|"                ${rundir}/input.geos
sed -i -e "s|{NESTED_SUM}|F|"                ${rundir}/input.geos
sed -i -e "s|{BUFFER_ZONE}|0  0  0  0|"      ${rundir}/input.geos
sed -i -e "s|{DATA_ROOT}|${GC_DATA_ROOT}|"   ${rundir}/HEMCO_Config.rc
sed -i -e "s|{GRID_DIR}|${grid_res_long}|"   ${rundir}/HEMCO_Config.rc
sed -i -e "s|{MET_DIR}|${met_dir}|"          ${rundir}/HEMCO_Config.rc
sed -i -e "s|{NATIVE_RES}|${met_native}|"    ${rundir}/HEMCO_Config.rc
sed -i -e "s|{LATRES}|${met_latres}|"        ${rundir}/HEMCO_Config.rc
sed -i -e "s|{LONRES}|${met_lonres}|"        ${rundir}/HEMCO_Config.rc
sed -i -e "s|{DUST_SF}|${dust_sf}|"          ${rundir}/HEMCO_Config.rc

# Special handling for start/end date based on simulation so that
# start year/month/day matches default initial restart file.
if [[ ${sim_type} =~ "chem" ]]; then
    startdate="20190701"
    enddate="20190801"
elif [[ ${sim_name} = "benchmark" ]]; then
    startdate="20190701"
    enddate="20190801"
else
    startdate="20190101"
    enddate="20190201"
fi
starttime="000000"
endtime="000000"
sed -i -e "s|{DATE1}|${startdate}|"     ${rundir}/input.geos
sed -i -e "s|{DATE2}|${enddate}|"       ${rundir}/input.geos
sed -i -e "s|{TIME1}|${starttime}|"     ${rundir}/input.geos
sed -i -e "s|{TIME2}|${endtime}|"       ${rundir}/input.geos

#-----------------------------------------------------------------
# Set permissions
#-----------------------------------------------------------------
chmod 755 ${rundir}/setCodeDir.sh
chmod 755 ${rundir}/cleanRundir.sh
chmod 755 ${rundir}/archiveRun.sh
chmod 755 ${rundir}/runScriptSamples/*

#----------------------------------------------------------------------
# Archive GCHP repository version in run directory file rundir.version
#----------------------------------------------------------------------
version_log=${rundir}/rundir.version
echo "This run directory was created with GEOS-Chem/run/GCClassic/createRunDir.sh." > ${version_log}
echo " " >> ${version_log}
echo "GEOS-Chem repository version information:" >> ${version_log}
cd ${gcdir}
remote_url=$(git config --get remote.origin.url)
code_branch=$(git rev-parse --abbrev-ref HEAD)
last_commit=$(git log -n 1 --pretty=format:"%s")
commit_date=$(git log -n 1 --pretty=format:"%cd")
commit_user=$(git log -n 1 --pretty=format:"%cn")
commit_hash=$(git log -n 1 --pretty=format:"%h")
cd ${curdir}
printf "\n  Remote URL: ${remote_url}" >> ${version_log}
printf "\n  Branch: ${code_branch}"    >> ${version_log}
printf "\n  Commit: ${last_commit}"    >> ${version_log}
printf "\n  Date: ${commit_date}"      >> ${version_log}
printf "\n  User: ${commit_user}"      >> ${version_log}
printf "\n  Hash: ${commit_hash}"      >> ${version_log}

#-----------------------------------------------------------------
# Ask user whether to track run directory changes with git
#-----------------------------------------------------------------
printf "\nDo you want to track run directory changes with git? (y/n)\n"
valid_response=0
while [ "$valid_response" -eq 0 ]
do
    read enable_git
    if [[ ${enable_git} = "y" ]]; then
	cd ${rundir}
	printf "\n\nChanges to the following run directory files are tracked by git:\n\n" >> ${version_log}
	git init
	git add *.rc *.sh *.yml *.run input.geos getRunInfo
	git add runScriptSamples/* README .gitignore
	git add 
	printf " " >> ${version_log}
	git commit -m "Initial run directory" >> ${version_log}
	cd ${curdir}
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
