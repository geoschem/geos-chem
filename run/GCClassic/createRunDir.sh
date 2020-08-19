#!/bin/bash

# createRunDir.sh: Create GEOS-Chem Classic run directory
#
# Optional argument: run directory name
#
# If optional run directory name argument is not passed then the user
# will be prompted to enter a name interactively, or choose to use the
# default name GC_{res}_{met}_{simulation}.
#
# Usage: ./createRunDir.sh [rundirname]
#
# Initial version: M. Sulprizio, 6/24/2020 (based off GCHPctm/createRunDir.sh)

srcrundir=$(pwd -P)
cd ${srcrundir}
cd ../..
gcdir=$(pwd -P)
cd ../../
wrapperdir=$(pwd -P)
cd ${srcrundir}

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
    printf "${thinline}Define path to ExtData"
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
	read extdata
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
printf "   2. Aerosols only\n"
printf "   3. CH4\n"
printf "   4. CO2\n"
printf "   5. Hg\n"
printf "   6. POPs\n"
printf "   7. Tagged CH4\n"
printf "   8. Tagged CO\n"
printf "   9. Tagged O3\n"
printf "  10. TransportTracers\n"

valid_sim=0
while [ "${valid_sim}" -eq 0 ]; do
    read sim_num
    valid_sim=1
    if [[ ${sim_num} = "1" ]]; then
	sim_name=fullchem
    elif [[ ${sim_num} = "2" ]]; then
	sim_name=aerosol
    elif [[ ${sim_num} = "3" ]]; then
	sim_name=CH4
    elif [[ ${sim_num} = "4" ]]; then
	sim_name=CO2
    elif [[ ${sim_num} = "5" ]]; then
	sim_name=Hg
    elif [[ ${sim_num} = "6" ]]; then
	sim_name=POPs
    elif [[ ${sim_num} = "7" ]]; then
	sim_name=tagCH4
    elif [[ ${sim_num} = "8" ]]; then
	sim_name=tagCO
    elif [[ ${sim_num} = "9" ]]; then
	sim_name=tagO3
    elif [[ ${sim_num} = "10" ]]; then
	sim_name=TransportTracers
    else
        valid_sim=0
	printf "Invalid simulation option. Try again.\n"
    fi
done

#-----------------------------------------------------------------
# Ask user to specify full-chemistry simulation options
#-----------------------------------------------------------------
if [[ ${sim_name} = "fullchem" ]]; then
    
    printf "${thinline}Choose chemistry domain:${thinline}"
    printf "  1. Troposphere + stratosphere (Recommended)\n"
    printf "  2. Troposphere only\n"
    valid_chemgrid=0
    while [ "${valid_chemgrid}" -eq 0 ]; do
	read chemgrid_num
	if [[ ${chemgrid_num} = "1" ]]; then
	    chemgrid="trop+strat"
	    valid_chemgrid=1
	elif [[ ${chemgrid_num} = "2" ]]; then
	    chemgrid="trop_only"
	    valid_chemgrid=1
	else
	  printf "Invalid chemistry domain option. Try again.\n"
	fi
    done
    
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
    sim_extra_option=none
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
	else
	    valid_sim_option=0
	    printf "Invalid simulation option. Try again.\n"
	fi
    done
    
#-----------------------------------------------------------------
# Ask user to specify POPs simulation options
#-----------------------------------------------------------------
elif [[ ${sim_name} = "POPs" ]]; then
    printf "${thinline}Choose POPs type:${thinline}"
    printf "  1. BaP\n"
    printf "  2. PHE\n"
    printf "  3. PYR\n"
    valid_pops=0
    while [ "${valid_pops}" -eq 0 ]; do
	read pops_num
	valid_pops=1
	if [[ ${pops_num} = "1" ]]; then
	    POP_SPC="BaP"
	    POP_XMW="252.31d-3"
	    POP_KOA="3.02d11"
	    POP_KBC="7.94d13"
	    POP_K_POPG_OH="50d-12"
	    POP_K_POPG_O3A="0d0"
	    POP_K_POPG_O3B="2.8d15"
	    POP_HSTAR="3.10d-5"
	    POP_DEL_H="-110d3"
	    POP_DEL_Hw="43d0"
	elif [[ ${pops_num} = "2" ]]; then
	    POP_SPC="PHE"
	    POP_XMW="178.23d-3"
	    POP_KOA="4.37d7"
	    POP_KBC="1.0d10"
	    POP_K_POPG_OH="2.7d-11"
	    POP_K_POPG_O3A="0d0"
	    POP_K_POPG_O3B="2.15d15"
	    POP_HSTAR="1.74d-3"
	    POP_DEL_H="-74d3"
	    POP_DEL_Hw="47d0"
	elif [[ ${pops_num} = "3" ]]; then
	    POP_SPC="PYR"
	    POP_XMW="202.25d-3"
	    POP_KOA="7.24d8"
	    POP_KBC="1.0d11"
	    POP_K_POPG_OH="50d-12"
	    POP_K_POPG_O3A="0d0"
	    POP_K_POPG_O3B="3.0d15"
	    POP_HSTAR="5.37d-4"
	    POP_DEL_H="-87d3"
	    POP_DEL_Hw="43d0"
	else
	    valid_pops=0
	    printf "Invalid POPs type. Try again.\n"
	fi
	sim_name="${sim_name}_${pops_spc}"
    done

fi

#-----------------------------------------------------------------
# Ask user to select meteorology source
#-----------------------------------------------------------------
printf "${thinline}Choose meteorology source:${thinline}"
printf "  1. MERRA2 (Recommended)\n"
printf "  2. GEOS-FP \n"
valid_met=0
while [ "${valid_met}" -eq 0 ]; do
    read met_num
    valid_met=1
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
	dust_sf='6.42e-5'
    else
	valid_met=0
	printf "Invalid meteorology option. Try again.\n"
    fi
done

#-----------------------------------------------------------------
# Ask user to select horizontal resolution
#-----------------------------------------------------------------
printf "${thinline}Choose horizontal resolution:${thinline}"
printf "  1. 4.0  x 5.0\n"
printf "  2. 2.0  x 2.5\n"
printf "  3. 0.5  x 0.625\n"
if [[ ${met_name} = "GEOSFP" ]]; then
    printf "  4. 0.25 x 0.3125\n"
fi

valid_res=0
while [ "${valid_res}" -eq 0 ]; do
    read res_num
    valid_res=1
    if [[ ${res_num} = "1" ]]; then
	grid_res='4x5'
	grid_res_long='4.0x5.0'
	grid_dir=$grid_res
    elif [[ ${res_num} = "2" ]]; then
	grid_res='2x25'
	grid_res_long='2.0x2.5'
	grid_dir='2x2.5'
    elif [[ ${res_num} = "3" ]]; then
	grid_res='05x0625'
	grid_res_long='0.5x0.625'
	grid_dir=$grid_res_long
    elif [[ ${res_num} = "4" ]]; then
	grid_res='025x03125'
	grid_res_long='0.25x0.3125'
	grid_dir=$grid_res_long
    else
	valid_res=0
	printf "Invalid horizontal resolution option. Try again.\n"
    fi
done

if [[ ${grid_res} = "05x0625" ]] || [[ ${grid_res} = "025x03125" ]]; then
    printf "${thinline}Choose horizontal grid domain:${thinline}"
    printf "  1. Global\n"
    printf "  2. Asia\n"
    printf "  3. Europe\n"
    printf "  4. North America\n"
    printf "  5. Custom\n"

    valid_domain=0
    while [ "${valid_domain}" -eq 0 ]; do
	read domain_num
	valid_domain=1
	if [[ ${domain_num} = "1" ]]; then
	    domain_name="global"
	    lon_range="-180.0 180.0"
	    lat_range=" -90.0  90.0"
	    half_polar="T"
	    nested_sim="F"
	    buffer_zone="0  0  0  0"
	else
	    domain_name="AS"
	    half_polar="F"
	    nested_sim="T"
	    buffer_zone="3  3  3  3"
	    if [[ ${domain_num} = "2" ]]; then
	        if [[ ${grid_res} = "05x0625" ]]; then 
	            lon_range=" 60.0 150.0"
		    lat_range="-11.0  55.0"
		elif [[ ${grid_res} = "025x03125" ]]; then 
	            lon_range=" 70.0 140.0"
		    lat_range=" 15.0  55.0"
		fi
	    elif [[ ${domain_num} = "3" ]]; then
		domain_name="EU"
	        if [[ ${grid_res} = "05x0625" ]]; then 
	            lon_range="-30.0 50.0"
		    lat_range=" 30.0 70.0"
		elif [[ ${grid_res} = "025x03125" ]]; then 
	            lon_range="-15.0  40.0"
		    lat_range=" 32.75 61.25"
		fi
	    elif [[ ${domain_num} = "4" ]]; then
		domain_name="NA"
	        if [[ ${grid_res} = "05x0625" ]]; then 
	            lon_range="-140.0 -40.0"
		    lat_range="  10.0  70.0"
		elif [[ ${grid_res} = "025x03125" ]]; then 
	            lon_range="-130.0  -60.0"
		    lat_range="   9.75  60.0"
		fi
	    elif [[ ${domain_num} = "5" ]]; then
		domain_name="custom"
	        lon_range="MinLon MaxLon"
	        lat_range="MinLat MaxLat"
	        printf "\n  -- You will need to manually set longitude and latitude"
		printf "\n     bounds in the Grid Menu of input.geos.\n"
	    else
  		valid_domain=0
		printf "Invalid horizontal grid domain option. Try again.\n"
	    fi
        fi
    done
else
    lon_range="-180.0 180.0"
    lat_range=" -90.0  90.0"
    half_polar="T"
    nested_sim="F"
    buffer_zone="0  0  0  0"
fi

#-----------------------------------------------------------------
# Ask user to select vertical resolution
#-----------------------------------------------------------------
printf "${thinline}Choose number of levels:${thinline}"
printf "  1. 72 (native)\n"
printf "  2. 47 (reduced)\n"

valid_lev=0
while [ "${valid_lev}" -eq 0 ]; do
    read lev_num
    valid_lev=1
    if [[ ${lev_num} = "1" ]]; then
	grid_lev='72'
    elif [[ ${lev_num} = "2" ]]; then
	grid_lev='47'
    else
	valid_lev=0
	printf "Invalid vertical resolution option. Try again.\n"
    fi
done

#-----------------------------------------------------------------
# Ask user to define path where directoy will be created
#-----------------------------------------------------------------
printf "${thinline}Enter path where the run directory will be created:${thinline}"
valid_path=0
while [ "$valid_path" -eq 0 ]; do
    read rundir_path
    if [[ ${rundir_path} = "q" ]]; then
	printf "\nExiting.\n"
	exit 1
    elif [[ ! -d ${rundir_path} ]]; then
        printf "\nERROR: ${rundir_path} does not exist. Enter a new path or hit q to quit.\n"
    else
	valid_path=1
    fi
done

#-----------------------------------------------------------------
# Ask user to define run directoy name if not passed as argument
#-----------------------------------------------------------------
if [ -z "$1" ]; then
    printf "${thinline}Enter run directory name, or press return to use default:${thinline}"
    read rundir_name
    if [[ -z "${rundir_name}" ]]; then
	if [[ "${sim_extra_option}" == "none" ]]; then
	    rundir_name=GC_${grid_res}_${sim_name}
	else
	    rundir_name=GC_${grid_res}_${sim_name}_${sim_extra_option}
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
	printf "WARNING: ${rundir} already exists.\n"
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
#------------------------------------------------------------------
cp -r ./getRunInfo  ${rundir}
cp -r ${gcdir}/run/shared/download_data.py  ${rundir}
cp -r ./runScriptSamples  ${rundir}
cp ./archiveRun.sh   ${rundir}
cp ./cleanRundir.sh  ${rundir}
cp ./setCodeDir.sh  ${rundir}
cp ./README  ${rundir}
cp ./gitignore  ${rundir}/.gitignore
cp ./input.geos.templates/input.geos.${sim_name}  ${rundir}/input.geos
cp ./HISTORY.rc.templates/HISTORY.rc.${sim_name}  ${rundir}/HISTORY.rc
cp ./HEMCO_Config.rc.templates/HEMCO_Config.rc.${sim_name}  ${rundir}/HEMCO_Config.rc
if [[ ${sim_name} =~ "chem" ]]; then
    cp ./HEMCO_Diagn.rc.templates/HEMCO_Diagn.rc.${sim_name}  ${rundir}/HEMCO_Diagn.rc
    cp -r ${gcdir}/run/shared/metrics_fullchem.py  ${rundir}
    chmod 744 ${rundir}/metrics_fullchem.py
else
    cp ./HEMCO_Diagn.rc.templates/HEMCO_Diagn.rc.${sim_name}  ${rundir}/HEMCO_Diagn.rc
fi
mkdir ${rundir}/OutputDir

# Copy species database; append APM or TOMAS species if needed
cp -r ${gcdir}/run/shared/species_database.yml   ${rundir}
if [[ ${sim_extra_option} =~ "TOMAS" ]]; then
    cat ${gcdir}/run/shared/species_database_tomas.yml >> ${rundir}/species_database.yml
elif [[ ${sim_extra_option} =~ "APM" ]]; then
    cat ${gcdir}/run/shared/species_database_apm.yml >> ${rundir}/species_database.yml
fi

# If benchmark simulation, put run script in directory
if [ ${sim_extra_option} = "benchmark" ]; then
    cp ./runScriptSamples/geoschem.benchmark.run ${rundir}
    chmod 744 ${rundir}/geoschem.benchmark.run
fi

#--------------------------------------------------------------------
# Create symbolic link to code directory
#--------------------------------------------------------------------
ln -s ${wrapperdir} ${rundir}/CodeDir

#--------------------------------------------------------------------
# Create build directory
#--------------------------------------------------------------------
mkdir ${rundir}/build
printf "To build GEOS-Chem type:\n   cmake ../CodeDir\n   make -j\n   make install\n" >> ${rundir}/build/README

#-----------------------------------------------------------------
# Replace token strings in certain files
#-----------------------------------------------------------------
sed -i -e "s|{DATA_ROOT}|${GC_DATA_ROOT}|"   ${rundir}/input.geos
sed -i -e "s|{MET}|${met_name}|"             ${rundir}/input.geos
sed -i -e "s|{SIM}|${sim_name}|"             ${rundir}/input.geos
sed -i -e "s|{RES}|${grid_res_long}|"        ${rundir}/input.geos
sed -i -e "s|{NLEV}|${grid_lev}|"            ${rundir}/input.geos
sed -i -e "s|{LON_RANGE}|${lon_range}|"      ${rundir}/input.geos
sed -i -e "s|{LAT_RANGE}|${lat_range}|"      ${rundir}/input.geos
sed -i -e "s|{HALF_POLAR}|${half_polar}|"    ${rundir}/input.geos
sed -i -e "s|{NESTED_SIM}|${nested_sim}|"    ${rundir}/input.geos
sed -i -e "s|{BUFFER_ZONE}|${buffer_zone}|"  ${rundir}/input.geos
sed -i -e "s|{DATA_ROOT}|${GC_DATA_ROOT}|"   ${rundir}/HEMCO_Config.rc
sed -i -e "s|{GRID_DIR}|${grid_dir}|"        ${rundir}/HEMCO_Config.rc
sed -i -e "s|{MET_DIR}|${met_dir}|"          ${rundir}/HEMCO_Config.rc
sed -i -e "s|{NATIVE_RES}|${met_native}|"    ${rundir}/HEMCO_Config.rc
sed -i -e "s|{LATRES}|${met_latres}|"        ${rundir}/HEMCO_Config.rc
sed -i -e "s|{LONRES}|${met_lonres}|"        ${rundir}/HEMCO_Config.rc
sed -i -e "s|{DUST_SF}|${dust_sf}|"          ${rundir}/HEMCO_Config.rc

#-----------------------------------------------------------------
# Update config file default settings based on simulation selected
#-----------------------------------------------------------------

#### Define function to replace values in config files
replace_colon_sep_val() {
    KEY=$1
    VALUE=$2
    FILE=$3
    #printf '%-30s : %-20s %-20s\n' "${KEY//\\}" "${VALUE}" "${FILE}"

    # replace value in line starting with 'whitespace + key + whitespace + : +
    # whitespace + value' where whitespace is variable length including none
    sed -i -e "s|^\([\t ]*${KEY}[\t ]*:[\t ]*\).*|\1${VALUE}|" ${FILE}
}

#### Define function to remove line(s) in config files
remove_text() {
    VALUE=$1
    FILE=$2
    sed -i -e "/${VALUE}/d" ${FILE}
}

#-----------------------------
# Benchmark settings
#-----------------------------
if [[ ${sim_extra_option} = "benchmark" ]]; then
    replace_colon_sep_val "Use GC classic timers?"  T     ${rundir}/input.geos
    replace_colon_sep_val "--> OFFLINE_DUST"        false ${rundir}/HEMCO_Config.rc
    replace_colon_sep_val "--> OFFLINE_BIOGENICVOC" false ${rundir}/HEMCO_Config.rc
    replace_colon_sep_val "--> OFFLINE_SEASALT"     false ${rundir}/HEMCO_Config.rc
    replace_colon_sep_val "--> OFFLINE_SOILNOX"     false ${rundir}/HEMCO_Config.rc
    sed -i -e "s|DustDead               : off|DustDead               : on |" ${rundir}/HEMCO_Config.rc
    sed -i -e "s|SoilNOx                : off|SoilNOx                : on |" ${rundir}/HEMCO_Config.rc
    sed -i -e "s|SeaSalt                : off|SeaSalt                : on |" ${rundir}/HEMCO_Config.rc
    sed -i -e "s|NO     0      3 |NO     104    -1|" ${rundir}/HEMCO_Diagn.rc
    sed -i -e "s|0      3 |105    -1|"               ${rundir}/HEMCO_Diagn.rc
    sed -i -e "s|0      4 |108    -1|"               ${rundir}/HEMCO_Diagn.rc
    sed -i -e "s|#Inv|Inv|"                          ${rundir}/HEMCO_Diagn.rc
    sed -i -e "s|#'|'|"                              ${rundir}/HISTORY.rc
fi

#-----------------------------
# Complex SOA settings
#-----------------------------
if [[ ${sim_extra_option} = "benchmark"   ]] || \
   [[ ${sim_extra_option} =~ "complexSOA" ]] || \
   [[ ${sim_extra_option} = "APM" ]]; then

    # Turn on complex SOA option in input.geos
    replace_colon_sep_val "Online COMPLEX SOA" T ${rundir}/input.geos

    # Add complex SOA species in input.geos
    prev_line="Species name            : ALK4"
    new_line="\Species name            : ASOA1\n\
Species name            : ASOA2\n\
Species name            : ASOA3\n\
Species name            : ASOAN\n\
Species name            : ASOG1\n\
Species name            : ASOG2\n\
Species name            : ASOG3"
    sed -i -e "/${prev_line}/a ${new_line}" ${rundir}/input.geos
    
    prev_line="Species name            : TOLU"
    new_line="\Species name            : TSOA0\n\
Species name            : TSOA1\n\
Species name            : TSOA2\n\
Species name            : TSOA3\n\
Species name            : TSOG0\n\
Species name            : TSOG1\n\
Species name            : TSOG2\n\
Species name            : TSOG3"
    sed -i -e "/${prev_line}/a ${new_line}" ${rundir}/input.geos
fi

#-----------------------------
# Semivolatile POA settings
#-----------------------------
if [[ ${sim_extra_option} = "complexSOA_SVPOA" ]]; then
	
    # Turn on semivolatile POA option in input.geos
    replace_colon_sep_val "=> Semivolatile POA?" T ${rundir}/input.geos

    # Add semivolatile POA species in input.geos
    prev_line="Species name            : N2O5"
    new_line="\Species name            : NAP"
    sed -i -e "/${prev_line}/a ${new_line}" ${rundir}/input.geos
	
    line="Species name            : OCPI"
    remove_text $line ${rundir}/input.geos
    line="Species name            : OCPO"
    remove_text $line ${rundir}/input.geos
	
    prev_line="Species name            : OIO"
    new_line="\Species name            : OPOA1\n\
Species name            : OPOA2\n\
Species name            : OPOG1\n\
Species name            : OPOG2"
    sed -i -e "/${prev_line}/a ${new_line}" ${rundir}/input.geos
	
    prev_line="Species name            : PIP"
    new_line="\Species name            : POA1\n\
Species name            : POA2\n\
Species name            : POG1\n\
Species name            : POG2"
    sed -i -e "/${prev_line}/a ${new_line}" ${rundir}/input.geos
fi

#-----------------------------
# Acid uptake settings
#-----------------------------
if [[ ${sim_extra_option} = "aciduptake" ]]; then
    replace_colon_sep_val "DustAlk"          on ${rundir}/HEMCO_Config.rc
    replace_colon_sep_val "=> Acidic uptake" T  ${rundir}/input.geos

    # Add acid uptake species in input.geos
    prev_line="Species name            : DST4"
    new_line="\Species name            : DSTAL1\n\
Species name            : DSTAL2\n\
Species name            : DSTAL3\n\
Species name            : DSTAL4"
    sed -i -e "/${prev_line}/a ${new_line}" ${rundir}/input.geos
    
    prev_line="Species name            : NIT"
    new_line="\Species name            : NITD1\n\
Species name            : NITD2\n\
Species name            : NITD3\n\
Species name            : NITD4"
    sed -i -e "/${prev_line}/a ${new_line}" ${rundir}/input.geos
    
    prev_line="Species name            : SO4"
    new_line="\Species name            : SO4D1\n\
Species name            : SO4D2\n\
Species name            : SO4D3\n\
Species name            : SO4D4"
    sed -i -e "/${prev_line}/a ${new_line}" ${rundir}/input.geos
fi

#-----------------------------
# Marine POA settings
#-----------------------------
if [[ ${sim_extra_option} = "marinePOA" ]]; then
    replace_colon_sep_val "SeaSalt"                 on ${rundir}/HEMCO_Config.rc
    replace_colon_sep_val " => MARINE ORG AEROSOLS" T  ${rundir}/input.geos

    # Add marine POA species to input.geos
    prev_line"Species name            : MONITU"
    new_line="\Species name            : MOPI\n\
Species name            : MOPO"
    sed -i -e "/${prev_line}/a ${new_line}" ${rundir}/input.geos
fi

#-----------------------------
# RRTMG settings
#-----------------------------
if [[ ${sim_extra_option} = "RRTMG" ]]; then

    replace_colon_sep_val "Turn on RRTMG?"       T ${rundir}/input.geos
    replace_colon_sep_val "Calculate LW fluxes?" T ${rundir}/input.geos
    replace_colon_sep_val "Calculate SW fluxes?" T ${rundir}/input.geos
    replace_colon_sep_val "Clear-sky flux?"      T ${rundir}/input.geos
    replace_colon_sep_val "All-sky flux?"        T ${rundir}/input.geos
    replace_colon_sep_val "--> RRTMG"         true ${rundir}/HEMCO_Config.rc
    sed -i -e "s|##'RRTMG'|'RRTMG'|"                ${rundir}/HISTORY.rc
fi

#-----------------------------
# TOMAS settings
#-----------------------------
if [[ ${sim_extra_option} =~ "TOMAS" ]]; then
    replace_colon_sep_val "Tran/conv timestep [sec]" 1800 ${rundir}/input.geos
    replace_colon_sep_val "Chem/emis timestep [sec]" 3600 ${rundir}/input.geos
    replace_colon_sep_val "Use non-local PBL?"       F    ${rundir}/input.geos
    replace_colon_sep_val "Use linear. strat. chem?" F    ${rundir}/input.geos
    replace_colon_sep_val "=> Online O3 from model"  F    ${rundir}/input.geos
    replace_colon_sep_val "TOMAS_Jeagle"             on   ${rundir}HEMCO_Config.rc
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
Species name            : ECOB1\n\
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
    sed -i -e "/${prev_line}/a ${new_line}" ${rundir}/input.geos
    
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
	sed -i -e "/${prev_line}/a ${new_line}" ${rundir}/input.geos

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
	sed -i -e "/${prev_line}/a ${new_line}" ${rundir}/input.geos

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
	sed -i -e "/${prev_line}/a ${new_line}" ${rundir}/input.geos

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
	sed -i -e "/${prev_line}/a ${new_line}" ${rundir}/input.geos

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
	sed -i -e "/${prev_line}/a ${new_line}" ${rundir}/input.geos

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
	sed -i -e "/${prev_line}/a ${new_line}" ${rundir}/input.geos

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
	sed -i -e "/${prev_line}/a ${new_line}" ${rundir}/input.geos

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
	sed -i -e "/${prev_line}/a ${new_line}" ${rundir}/input.geos

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
	sed -i -e "/${prev_line}/a ${new_line}" ${rundir}/input.geos

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
    sed -i -e "/${prev_line}/a ${new_line}" ${rundir}/input.geos
    
fi
    
#-----------------------------------------------------------------
# Modify input files for troposphere-only chemistry grids
#-----------------------------------------------------------------
if [[ ${chemgrid} = "trop_only" ]]; then

    replace_colon_sep_val "=> Set init. strat. H2O"  F ${rundir}/input.geos
    replace_colon_sep_val "Settle strat. aerosols"   F ${rundir}/input.geos
    replace_colon_sep_val "Online PSC AEROSOLS"      F ${rundir}/input.geos
    replace_colon_sep_val "Perform PSC het. chem.?"  F ${rundir}/input.geos
    replace_colon_sep_val "Calc. strat. aero. OD?"   F ${rundir}/input.geos
    replace_colon_sep_val "Use UCX strat. chem?"     F ${rundir}/input.geos
    replace_colon_sep_val "Active strat. H2O?"       F ${rundir}/input.geos
    replace_colon_sep_val "--> STATE_PSC"        false ${rundir}/HEMCO_Config.rc
    replace_colon_sep_val "--> GMI_PROD_LOSS"    false ${rundir}/HEMCO_Config.rc
    replace_colon_sep_val "--> UCX_PROD_LOSS"     true ${rundir}/HEMCO_Config.rc
    sed -i -e "s|'Chem_StatePSC|#'Chem_StatePSC|"      ${rundir}/HISTORY.rc

fi

#-----------------------------------------------------------------
# Modify input files for nested-grid simulations
#-----------------------------------------------------------------
if [[ ${nested_sim} = "T" ]]; then
    replace_colon_sep_val "--> GC_BCs" true ${rundir}/HEMCO_Config.rc
    if [[ ${domain_name} = "NA" ]]; then
	replace_colon_sep_val "--> NEI2011_MONMEAN" false ${rundir}/HEMCO_Config.rc
	replace_colon_sep_val "--> NEI2011_HOURLY"  true  ${rundir}/HEMCO_Config.rc
    fi
fi

#-----------------------------------------------------------------
# Modify input files for POPs simulations
#-----------------------------------------------------------------
if [[ ${sim_name} =~ "POPs" ]]; then
    sed -i -e "s|{POPs_SPC}|${POP_SPC}|"               ${rundir}/input.geos
    sed -i -e "s|{POPs_XMW}|${POP_XMW}|"               ${rundir}/input.geos
    sed -i -e "s|{POPs_KOA}|${POP_KOA}|"               ${rundir}/input.geos
    sed -i -e "s|{POPs_KBC}|${POP_KBC}|"               ${rundir}/input.geos
    sed -i -e "s|{POPs_K_POPG_OH}|${POP_K_POPG_OH}|"   ${rundir}/input.geos
    sed -i -e "s|{POPs_K_POPG_O3A}|${POP_K_POPG_O3A}|" ${rundir}/input.geos
    sed -i -e "s|{POPs_K_POPG_O3B}|${POP_K_POPG_O3B}|" ${rundir}/input.geos
    sed -i -e "s|{POPs_HSTAR}|${POP_HSTAR}|"           ${rundir}/input.geos
    sed -i -e "s|{POPs_DEL_H}|${POP_DEL_H}|"           ${rundir}/input.geos
    sed -i -e "s|{POPs_DEL_Hw}|${POP_DEL_Hw}|"         ${rundir}/input.geos
    sed -i -e "s|{POPs_SPC}|${POP_SPC}|"               ${rundir}/HEMCO_Config.rc
    sed -i -e "s|{POPs_SPC}|${POP_SPC}|"               ${rundir}/HEMCO_Diagn.rc
fi

#-----------------------------------------------------------------
# Special handling for start/end date based on simulation so that
# start year/month/day matches default initial restart file.
#-----------------------------------------------------------------
if [[ ${sim_extra_option} = "benchmark" ]]; then
    startdate="20190701"
    enddate="20190801"
elif [[ ${sim_name} = "fullchem" ]]; then
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

printf "\n  -- This run directory has been set up for $startdate - $enddate."
printf "\n     You may modify these settings in input.geos.\n"

sed -i -e "s|{FREQUENCY}|00000100 000000|"   ${rundir}/HISTORY.rc
sed -i -e "s|{DURATION}|00000100 000000|"    ${rundir}/HISTORY.rc

printf "\n  -- The default frequency and duration of diagnostics is set to monthly."
printf "\n     You may modify these settings in HISTORY.rc and HEMCO_Config.rc.\n"

#--------------------------------------------------------------------
# Copy sample restart file
#--------------------------------------------------------------------
if [[ ${sim_name} = "fullchem" ]]; then
    # Use restart file saved out from latest 1-year benchmark
    sample_rst=${GC_DATA_ROOT}/GEOSCHEM_RESTARTS/GC_12.9.0/initial_GEOSChem_rst.4x5_benchmark.nc
else
    sample_rst=${GC_DATA_ROOT}/GEOSCHEM_RESTARTS/v2018-11/initial_GEOSChem_rst.${grid_res}_${sim_name}.nc
fi
if [[ -f ${sample_rst} ]]; then
    cp ${sample_rst} ${rundir}/GEOSChem.Restart.${startdate}_0000z.nc4
else
    printf "\n  -- No sample restart provided for this simulation."
    printf "\n     You will need to provide an initial restart file or disable"
    printf "\n     GC_RESTARTS in HEMCO_Config.rc to initialize your simulation"
    printf "\n     with default background species concentrations.\n"
fi

#-----------------------------------------------------------------
# Set permissions
#-----------------------------------------------------------------
chmod 755 ${rundir}/setCodeDir.sh
chmod 755 ${rundir}/cleanRundir.sh
chmod 755 ${rundir}/archiveRun.sh
chmod 755 ${rundir}/runScriptSamples/*

#----------------------------------------------------------------------
# Archive repository version in run directory file rundir.version
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
	git add *.rc *.sh *.yml *.run *.py input.geos getRunInfo
	git add runScriptSamples/* README .gitignore
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
