#!/bin/bash

# createRunDir.sh: Create HEMCO standalone run directory
#
# Optional argument: run directory name
#
# If optional run directory name argument is not passed then the user
# will be prompted to enter a name interactively, or choose to use the
# default name gchp_{simulation}/
#
# Usage: ./createRunDir.sh [rundirname]
#
# Initial version: E. Lundgren,10/5/2018

curdir=$(pwd)
cd ..
hemcodir=$(pwd)
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
    if [[ ! -d ${GFTL} ]]; then
	printf "\nWarning: Default Goddard Fortran Template Library (gFTL) does not exist!"
        printf "\nSet new path below or manually edit ${HOME}/.geoschem/config.\n"
    fi
else
    printf "\nDefine paths to ExtData and the Goddard Fortran Template Library (gFTL)."
    printf "\nThese will be stored in ${HOME}/.geoschem/config for future automatic use.\n"
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
# Ask user to select meteorology source
#-----------------------------------------------------------------
printf "\nChoose meteorology source:\n"
printf "  1. GEOS-FP\n"
printf "  2. MERRA2\n"
valid_met=0
while [ "${valid_met}" -eq 0 ]
do
    read met_num
    if [[ ${met_num} = "1" ]]; then
	met_name='GEOSFP'
	met_dir='GEOS_FP'
	met_native='0.25x0.3125'
	met_latres='025'
	met_lonres='03125'
	met_resolution='${met_latres}x${met_lonres}'
	met_extension='nc'
	met_cn_year='2011'
	pressure_unit='hPa'
	pressure_scale='1.0 '
	valid_met=1
	dust_sf='6.42e-5'
    elif [[ ${met_num} = "2" ]]; then
	met_name='MERRA2'
	met_dir='MERRA2'
	met_native='0.5x0.625'
	met_latres='05'
	met_lonres='0625'
	met_resolution='${met_latres}x${met_lonres}'
	met_extension='nc4'
	met_cn_year='2015'
	pressure_unit='Pa '
	pressure_scale='0.01'
	valid_met=1
	dust_sf='3.86e-4'
    else
	printf "Invalid meteorology option. Try again.\n"
    fi
done

#-----------------------------------------------------------------
# Ask user to select meteorology source
#-----------------------------------------------------------------
printf "\nChoose horizontal resolution:\n"
printf "  1. 4x5\n"
printf "  2. 2x2.5\n"
printf "  3. 0.5x0.625\n"
printf "  4. 0.25x0.3125\n"
printf "  5. Custom\n"
valid_res=0
while [ "${valid_res}" -eq 0 ]
do
    read res_num
    if [[ ${res_num} = "1" ]]; then
	grid_name='4x5'
	grid_res='4x5'
	grid_file='HEMCO_sa_Grid.4x5.rc'
	valid_res=1
    elif [[ ${res_num} = "2" ]]; then
	grid_name='2x25'
	grid_res='2x2.5'
	grid_file='HEMCO_sa_Grid.2x25.rc'
	valid_res=1
    elif [[ ${res_num} = "3" ]]; then
	grid_name='05x0625'
	grid_res='0.5x0.625'
	grid_file='HEMCO_sa_Grid.05x0625.rc'
	valid_res=1
    elif [[ ${res_num} = "4" ]]; then
	grid_name='025x03125'
	grid_res='0.25x0.3125'
	grid_file='HEMCO_sa_Grid.025x03125.rc'
	valid_res=1
    elif [[ ${res_num} = "5" ]]; then
	printf "You will need to provide your own HEMCO_sa_Grid.rc file.\n"
	printf "See the HEMCO standalone guide for more information:\n"
	printf "http://wiki.seas.harvard.edu/geos-chem/index.php/HEMCO_standalone\n"
	valid_res=1
    else
	printf "Invalid resolution option. Try again.\n"
    fi
done

#-----------------------------------------------------------------
# Ask user to provide path to HEMCO_Config.template file
#-----------------------------------------------------------------
printf "\nEnter path to the HEMCO_Config.rc file with your emissions settings.\n\n"
printf "NOTE: This may be a HEMCO_Config.rc file from a GEOS-Chem run directory\n"
printf "or a HEMCO_Config.template file from the GEOS-Chem source code repository.\n"
valid_path=0
while [ "$valid_path" -eq 0 ]
do
    read hco_config_path
    if [[ ${hco_config_path} = "q" ]]; then
	printf "\nExiting.\n"
	exit 1
    elif [[ ! -f ${hco_config_path} ]]; then
        printf "\nError: ${hco_config_path} does not exist. Enter a new path or hit q to quit.\n"
    else
	valid_path=1
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
	rundir_name=hemco_${grid_name}_"${met_name,,}"
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
cp -r ./OutputDir ${rundir}
cp ${hco_config_path}          ${rundir}/HEMCO_Config.rc
cp ./HEMCO_sa_Config.template  ${rundir}/HEMCO_sa_Config.rc
cp ./HEMCO_sa_Time.rc          ${rundir}
cp ./HEMCO_sa_Spec.rc          ${rundir}
cp ./${grid_file}              ${rundir}
cp ./HEMCO_Diagn.rc            ${rundir}
cp ./runHEMCO.sh               ${rundir}
cp ./README                    ${rundir}
mkdir ${rundir}/build
printf "To build HEMCO type:\n   cmake ../CodeDir\n   make -j\n   make install\n" >> ${rundir}/build/README

#--------------------------------------------------------------------
# Create symbolic links to data directories, restart files, and code
#--------------------------------------------------------------------
ln -s ${hemcodir}                            ${rundir}/CodeDir

#-----------------------------------------------------------------
# Replace token strings in certain files
#-----------------------------------------------------------------
# HEMCO_sa_Config.rc
sed -i -e "s|{DATA_ROOT}|${GC_DATA_ROOT}|"   ${rundir}/HEMCO_sa_Config.rc
sed -i -e "s|{GRID_FILE}|${grid_file}|"      ${rundir}/HEMCO_sa_Config.rc
sed -i -e "s|{MET_NAME}|${met_name}|"        ${rundir}/HEMCO_sa_Config.rc
sed -i -e "s|{MET_RES}|${grid_res}|"         ${rundir}/HEMCO_sa_Config.rc

# HEMCO_Config.rc (copied from GEOS-Chem)
sed -i -e "s|{DATA_ROOT}|${GC_DATA_ROOT}|"   ${rundir}/HEMCO_Config.rc
sed -i -e "s|{GRID_DIR}|${grid_res}|"        ${rundir}/HEMCO_Config.rc
sed -i -e "s|{MET_DIR}|${met_dir}|"          ${rundir}/HEMCO_Config.rc
sed -i -e "s|{VERBOSE}|0|"                   ${rundir}/HEMCO_Config.rc
sed -i -e "s|{WARNINGS}|1|"                  ${rundir}/HEMCO_Config.rc
sed -i -e "s|{NATIVE_RES}|${met_native}|"    ${rundir}/HEMCO_Config.rc
sed -i -e "s|{LATRES}|${met_latres}|"        ${rundir}/HEMCO_Config.rc
sed -i -e "s|{LONRES}|${met_lonres}|"        ${rundir}/HEMCO_Config.rc
sed -i -e "s|{DUST_SF}|${dust_sf}|"          ${rundir}/HEMCO_Config.rc

#----------------------------------------------------------------------
# Archive GCHP repository version in run directory file rundir.version
#----------------------------------------------------------------------
version_log=${rundir}/rundir.version
echo "This run directory was created with HEMCO/run/createRunDir.sh." > ${version_log}
echo " " >> ${version_log}
echo "HEMCO repository version information:" >> ${version_log}
cd ${hemcodir}
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
	git add *.rc *.sh
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
