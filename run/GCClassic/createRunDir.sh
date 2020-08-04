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

curdir=$(pwd)
cd ../../../../
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
    while [ "$valid_path" -eq 0 ]; do
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
    if [[ ${sim_num} = "1" ]]; then
	sim_type=fullchem
	valid_sim=1
    elif [[ ${sim_num} = "2" ]]; then
	sim_type=aerosol
	valid_sim=1
    elif [[ ${sim_num} = "3" ]]; then
	sim_type=CH4
	valid_sim=1
    elif [[ ${sim_num} = "4" ]]; then
	sim_type=CO2
	valid_sim=1
    elif [[ ${sim_num} = "5" ]]; then
	sim_type=Hg
	valid_sim=1
    elif [[ ${sim_num} = "6" ]]; then
	sim_type=POPs
	valid_sim=1
    elif [[ ${sim_num} = "7" ]]; then
	sim_type=tagCH4
	valid_sim=1
    elif [[ ${sim_num} = "8" ]]; then
	sim_type=tagCO
	valid_sim=1
    elif [[ ${sim_num} = "9" ]]; then
	sim_type=tagO3
	valid_sim=1
    elif [[ ${sim_num} = "10" ]]; then
	sim_type=TransportTracers
	valid_sim=1
    else
	printf "Invalid simulation option. Try again.\n"
    fi
done

#-----------------------------------------------------------------
# Ask user to specify full-chemistry simulation options
#-----------------------------------------------------------------
if [[ ${sim_type} = "fullchem" ]]; then
    
    printf "\nChoose chemistry domain:\n"
    printf "  1. Troposphere + stratosphere (recommended)\n"
    printf "  2. Troposphere only\n"
    valid_chemgrid=0
    while [ "${valid_chemgrid}" -eq 0 ]; do
	read chemgrid_num
	if [[ ${chemgrid_num} = "1" ]]; then
	    chemgrid="troposphere"
	    valid_chemgrid=1
	elif [[ ${chemgrid_num} = "2" ]]; then
	    chemgrid="tropstrat"
	    valid_chemgrid=1
	else
	  printf "Invalid chemistry domain option. Try again.\n"
	fi
    done
    
    printf "\nChoose additional simulation option:\n"
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
	if [[ ${sim_option} = "1" ]]; then
	    sim_name="standard"
	    valid_sim_option=1
	elif [[ ${sim_option} = "2" ]]; then
	    sim_name="benchmark"
	    valid_sim_option=1
	elif [[ ${sim_option} = "3" ]]; then
	    printf "\nChoose complex SOA option:\n"
	    printf "  1. Complex SOA\n"
	    printf "  2. Complex SOA with semivolatile POA\n"
	    valid_soa=0
	    while [ "${valid_soa}" -eq 0 ]; do
		read soa_option
		if [[ ${soa_option} = "1" ]]; then
		    sim_name="complexSOA"
		    valid_soa=1
		elif [[ ${soa_option} = "2" ]]; then
		    sim_name="complexSOA_SVPOA"
		    valid_soa=1
		fi
	    done
	    valid_sim_option=1
	elif [[ ${sim_option} = "4" ]]; then
	   sim_name="marinePOA"
	   valid_sim_option=1
	elif [[ ${sim_option} = "5" ]]; then
	   sim_name="aciduptake"
	   valid_sim_option=1
	elif [[ ${sim_option} = "6" ]]; then
	    printf "\nChoose TOMAS option:\n"
	    printf "  1. TOMAS15\n"
	    printf "  1. TOMAS40\n"
	    valid_tomas=0
	    while [ "${valid_tomas}" -eq 0 ]; do
		read tomas_option
		if [[ ${tomas_option} = "1" ]]; then
		    sim_name="TOMAS15"
		    valid_tomas=1
		elif [[ ${tomas_option} = "2" ]]; then
		    sim_name="TOMAS40"
		    valid_tomas=1
		fi
	    done
	    valid_sim_option=1
	    valid_sim_option=1
	elif [[ ${sim_option} = "7" ]]; then
	    sim_name="APM"
	    valid_sim_option=1
	elif [[ ${sim_option} = "8" ]]; then
	    sim_name="RRTMG"
	    valid_sim_option=1
	else
	    printf "Invalid simulation option. Try again.\n"
	fi
    done
    
#-----------------------------------------------------------------
# Ask user to specify POPs simulation options
#-----------------------------------------------------------------
elif [[ ${sim_type} = "POPs" ]]; then
    printf "\nChoose POPs type:\n"
    printf "  1. BaP\n"
    printf "  2. PHE\n"
    printf "  3. PYR\n"
    valid_pops=0
    while [ "${valid_pops}" -eq 0 ]; do
	read pops_num
	if [[ ${pops_num} = "1" ]]; then
	    sim_name="POPs_BaP"
	    valid_pops=1
	elif [[ ${pops_num} = "2" ]]; then
	    sim_name="POPs_PHE"
	    valid_pops=1
	elif [[ ${pops_num} = "3" ]]; then
	    sim_name="POPs_PYR"
	    valid_pops=1
	else
	    printf "Invalid POPs type. Try again.\n"
	fi
    done

else
    sim_name=${sim_type}
fi

#-----------------------------------------------------------------
# Ask user to select meteorology source
#-----------------------------------------------------------------
printf "\nChoose meteorology source:\n"
printf "  1. MERRA2 (Recommended)\n"
printf "  2. GEOS-FP \n"
valid_met=0
while [ "${valid_met}" -eq 0 ]; do
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
while [ "${valid_res}" -eq 0 ]; do
    read res_num
    if [[ ${res_num} = "1" ]]; then
	grid_res='4x5'
	grid_res_long='4.0x5.0'
	grid_dir=$grid_res
	valid_res=1
    elif [[ ${res_num} = "2" ]]; then
	grid_res='2x25'
	grid_res_long='2.0x2.5'
	grid_dir='2x2.5'
	valid_res=1
    elif [[ ${res_num} = "3" ]]; then
	grid_res='05x0625'
	grid_res_long='0.5x0.625'
	grid_dir=$grid_res_long
	valid_res=1
    elif [[ ${res_num} = "4" ]]; then
	grid_res='025x03125'
	grid_res_long='0.25x0.3125'
	grid_dir=$grid_res_long
	valid_res=1
    else
	printf "Invalid horizontal resolution option. Try again.\n"
    fi
done

if [[ ${grid_res} = "05x0625" ]] || [[ ${grid_res} = "025x03125" ]]; then
    printf "\nChoose horizontal grid domain:\n"
    printf "  1. Global\n"
    printf "  2. Asia\n"
    printf "  3. Europe\n"
    printf "  4. North America\n"
    printf "  5. Custom\n"

    valid_domain=0
    while [ "${valid_domain}" -eq 0 ]; do
	read domain_num
	if [[ ${domain_num} = "1" ]]; then
	    lon_range="-180.0 180.0"
	    lat_range=" -90.0  90.0"
	    half_polar="T"
	    nested_sim="F"
	    buffer_zone="0  0  0  0"
	    valid_domain=1
	else
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
  		valid_domain=1
	    elif [[ ${domain_num} = "3" ]]; then
	        if [[ ${grid_res} = "05x0625" ]]; then 
	            lon_range="-30.0 50.0"
		    lat_range=" 30.0 70.0"
		elif [[ ${grid_res} = "025x03125" ]]; then 
	            lon_range="-15.0  40.0"
		    lat_range=" 32.75 61.25"
		fi
  		valid_domain=1
	    elif [[ ${domain_num} = "4" ]]; then
	        if [[ ${grid_res} = "05x0625" ]]; then 
	            lon_range="-140.0 -40.0"
		    lat_range="  10.0  70.0"
		elif [[ ${grid_res} = "025x03125" ]]; then 
	            lon_range="-130.0  -60.0"
		    lat_range="   9.75  60.0"
		fi
  		valid_domain=1
	    elif [[ ${domain_num} = "5" ]]; then
	        lon_range="MinLon MaxLon"
	        lat_range="MinLat MaxLat"
  		valid_domain=1
	        printf "\n NOTE: You will need to manually set longitude and latitude bounds in input.geos.\n"
	    else
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
printf "\nChoose number of levels:\n"
printf "  1. 72 (native)\n"
printf "  2. 47 (reduced)\n"

valid_lev=0
while [ "${valid_lev}" -eq 0 ]; do
    read lev_num
    if [[ ${lev_num} = "1" ]]; then
	grid_lev='72'
	valid_lev=1
    elif [[ ${lev_num} = "2" ]]; then
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
while [ "$valid_path" -eq 0 ]; do
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
	printf " Using default directory name ${rundir_name}\n"
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
cp ./input.geos.templates/input.geos.${sim_type}            ${rundir}/input.geos
cp ./HISTORY.rc.templates/HISTORY.rc.${sim_type}            ${rundir}/HISTORY.rc
cp ./HEMCO_Config.rc.templates/HEMCO_Config.rc.${sim_type}  ${rundir}/HEMCO_Config.rc
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
# Create symbolic link to code directory
#--------------------------------------------------------------------
ln -s ${gcdir} ${rundir}/CodeDir

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
sed -i -e "s|{FREQUENCY}|00000100 000000|"   ${rundir}/HISTORY.rc
sed -i -e "s|{DURATION}|00000100 000000|"    ${rundir}/HISTORY.rc

# Make changes for specific simulations
if [[ ${sim_name} = "benchmark" ]]; then
    line1="--> OFFLINE_DUST           :       true"
    line2="--> OFFLINE_DUST           :       false"
    sed -i -e "s|${line1}|${line2}|" ${rundir}/HEMCO_Config.rc
    
    line1="--> OFFLINE_BIOGENICVOC    :       true"
    line2="--> OFFLINE_BIOGENICVOC    :       false"
    sed -i -e "s|${line1}|${line2}|" ${rundir}/HEMCO_Config.rc

    line1="--> OFFLINE_SEASALT        :       true"
    line2="--> OFFLINE_SEASALT        :       false"
    sed -i -e "s|${line1}|${line2}|" ${rundir}/HEMCO_Config.rc

    line1="--> OFFLINE_SOILNOX        :       true"
    line2="--> OFFLINE_SOILNOX        :       false"
    sed -i -e "s|${line1}|${line2}|" ${rundir}/HEMCO_Config.rc

    line1="SoilNOx                : off"
    line2="SoilNOx                : on "
    sed -i -e "s|${line1}|${line2}|" ${rundir}/HEMCO_Config.rc

    line1="DustDead               : off"
    line2="DustDead               : on "
    sed -i -e "s|${line1}|${line2}|" ${rundir}/HEMCO_Config.rc

    line1="SeaSalt                : off"
    line2="SeaSalt                : on "
    sed -i -e "s|${line1}|${line2}|" ${rundir}/HEMCO_Config.rc

    line1="Use GC classic timers?  : F"
    line2="Use GC classic timers?  : T"
    sed -i -e "s|${line1}|${line2}|" ${rundir}/input.geos
elif [[ ${sim_name} = "aciduptake" ]]; then
    line1="DustAlk                : off"
    line2="DustAlk                : on "
elif [[ ${sim_name} = "marinePOA" ]]; then
    line1="SeaSalt                : off"
    line2="SeaSalt                : on "
    sed -i -e "s|${line1}|${line2}|" ${rundir}/HEMCO_Config.rc
fi

# Add species for special simulation options
if [[ ${sim_name} = "benchmark" ]] || [[ ${sim_name} =~ "complexSOA" ]]; then
    line1="Online COMPLEX SOA      : F"
    line2="Online COMPLEX SOA      : T"
    sed -i -e "s|${line1}|${line2}|" ${rundir}/input.geos
    
    line1="Species name            : ALK4"
    line2="\Species name            : ASOA1\n\
Species name            : ASOA2\n\
Species name            : ASOA3\n\
Species name            : ASOAN\n\
Species name            : ASOG1\n\
Species name            : ASOG2\n\
Species name            : ASOG3"
    sed -i -e "/${line1}/a ${line2}" ${rundir}/input.geos

    line1="Species name            : TOLU"
    line2="\Species name            : TSOA0\n\
Species name            : TSOA1\n\
Species name            : TSOA2\n\
Species name            : TSOA3\n\
Species name            : TSOG0\n\
Species name            : TSOG1\n\
Species name            : TSOG2\n\
Species name            : TSOG3"
    sed -i -e "/${line1}/a ${line2}" ${rundir}/input.geos

    if [[ ${sim_name} = "complexSOA_SVPOA" ]]; then
	line1="=> Semivolatile POA?   : F"
	line2="=> Semivolatile POA?   : T"
	sed -i -e "s|${line1}|${line2}|" ${rundir}/input.geos

	line1="Species name            : N2O5"
	line2="Species name            : NAP"
	sed -i -e "/${line1}/a ${line2}" ${rundir}/input.geos

	line1="Species name            : OCPI"
	line2="Species name            : OCPO"
	sed -i -e "/{$line1}/d"  ${rundir}/input.geos
	sed -i -e "/{$line2}/d"  ${rundir}/input.geos

	line1="Species name            : OIO"
	line2="\Species name            : OPOA1\n\
Species name            : OPOA2\n\
Species name            : OPOG1\n\
Species name            : OPOG2"
	sed -i -e "/${line1}/a ${line2}" ${rundir}/input.geos

	line1="Species name            : PIP"
	line2="\Species name            : POA1\n\
Species name            : POA2\n\
Species name            : POG1\n\
Species name            : POG2"
	sed -i -e "/${line1}/a ${line2}" ${rundir}/input.geos
    fi

elif [[ ${sim_name} = "aciduptake" ]]; then

    line1="=> Acidic uptake ?     : F"
    line2="=> Acidic uptake ?     : T"
    sed -i -e "s|${line1}|${line2}|" ${rundir}/input.geos
    
    line1="Species name            : DST4"
    line2="\Species name            : DSTAL1\n\
Species name            : DSTAL2\n\
Species name            : DSTAL3\n\
Species name            : DSTAL4"
    sed -i -e "/${line1}/a ${line2}" ${rundir}/input.geos

    line1="Species name            : NIT"
    line2="\Species name            : NITD1\n\
Species name            : NITD2\n\
Species name            : NITD3\n\
Species name            : NITD4"
    sed -i -e "/${line1}/a ${line2}" ${rundir}/input.geos

        line1="Species name            : SO4"
    line2="\Species name            : SO4D1\n\
Species name            : SO4D2\n\
Species name            : SO4D3\n\
Species name            : SO4D4"
    sed -i -e "/${line1}/a ${line2}" ${rundir}/input.geos

elif [[ ${sim_name} = "marinePOA" ]]; then

    line1=" => MARINE ORG AEROSOLS : F"
    line2=" => MARINE ORG AEROSOLS : T"
    sed -i -e "s|${line1}|${line2}|" ${rundir}/input.geos
    
    line1="Species name            : MONITU"
    line2="\Species name            : MOPI\n\
Species name            : MOPO"
    sed -i -e "/${line1}/a ${line2}" ${rundir}/input.geos

fi

# Save additional diagnostics for benchmark simulations
if [[ ${sim_name} = "benchmark" ]]; then
    line1="#Inv"
    line2="Inv"
    sed -i -e "s|${line1}|${line2}|" ${rundir}/HEMCO_Diagn.rc

    line1="#'"
    line2="'"
    sed -i -e "s|${line1}|${line2}|" ${rundir}/HISTORY.rc
fi


# Special handling for start/end date based on simulation so that
# start year/month/day matches default initial restart file.
if [[ ${sim_name} = "benchmark" ]]; then
    startdate="20190701"
    enddate="20190801"
elif [[ ${sim_type} =~ "chem" ]]; then
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

#--------------------------------------------------------------------
# Copy sample restart file
#--------------------------------------------------------------------
sample_rst=${GC_DATA_ROOT}/GEOSCHEM_RESTARTS/v2018-11/initial_GEOSChem_rst.${grid_res}_${sim_name}.nc
if [[ -f ${sample_rst} ]]; then
    cp ${sample_rst} ${rundir}/GEOSChem.Restart.${startdate}_0000z.nc4
else
    printf "\n NOTE: No sample restart provided for this simulation. You will need to provide an initial restart file or disable GC_RESTARTS in HEMCO_Config.rc to initialize your simulation with default background species concentrations.\n"
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
