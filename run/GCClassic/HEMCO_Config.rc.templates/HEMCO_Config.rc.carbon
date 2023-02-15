#------------------------------------------------------------------------------
#                  Harmonized Emissions Component (HEMCO)                     !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: HEMCO_Config.rc
#
# !DESCRIPTION: Contains configuration information for HEMCO. Define the
#  emissions inventories and corresponding file paths here. Entire
#  configuration files can be inserted into this configuration file with
#  an '>>>include' statement, e.g. '>>>include HEMCO\_Config\_test.rc'
#  The settings of include-files will be ignored.
#\\
#\\
# !REMARKS:
#  This file has been customized for the carbon simulation.
#  See The HEMCO User's Guide for file details:
#    http://wiki.geos-chem.org/The_HEMCO_User%27s_Guide
#
# !REVISION HISTORY:
#  See https://github.com/geoschem/geos-chem for complete history
#EOP
#------------------------------------------------------------------------------
#BOC
###############################################################################
### BEGIN SECTION SETTINGS
###############################################################################

ROOT:                        ${RUNDIR_DATA_ROOT}/HEMCO
METDIR:                      ${RUNDIR_MET_DIR}
GCAP2SCENARIO:               ${RUNDIR_GCAP2_SCENARIO}
GCAPVERTRES:                 ${RUNDIR_GCAP2_VERTRES}
Logfile:                     *
DiagnFile:                   HEMCO_Diagn.rc
DiagnPrefix:                 ./OutputDir/HEMCO_diagnostics
DiagnFreq:                   Monthly
Wildcard:                    *
Separator:                   /
Unit tolerance:              1
Negative values:             2
Only unitless scale factors: false
Verbose:                     0
Warnings:                    1
Mask fractions:              false

### END SECTION SETTINGS ###

###############################################################################
### BEGIN SECTION EXTENSION SWITCHES
###############################################################################
# ExtNr ExtName                on/off  Species   Years avail.
0       Base                   : on    *
# ----- MAIN SWITCHES ---------------------------------------------------------
    --> EMISSIONS              :       true
    --> METEOROLOGY            :       true      ${RUNDIR_MET_AVAIL}
    --> CHEMISTRY_INPUT        :       true
    --> USE_CH4_DATA           :       true
    --> USE_CO_DATA            :       true
    --> USE_CO2_DATA           :       true
    --> USE_OCS_DATA           :       true     # 2012-2016
# ----- RESTART FIELDS --------------------------------------------------------
    --> GC_RESTART             :       true
    --> HEMCO_RESTART          :       true
# ----- NESTED GRID FIELDS ----------------------------------------------------
    --> GC_BCs                 :       ${RUNDIR_USE_BCs}
# ----- CH4 INVENTORIES AND DATA ----------------------------------------------
# ..... Regional Inventories .........
    --> GEPA                   :       true     # 2012
    --> Scarpelli_Canada       :       true     # 2018
    --> Scarpelli_Mexico       :       true     # 2015
# ..... Global Inventories ...........
    --> GFEIv2                 :       true     # 2019
    --> EDGARv6                :       true     # 2000-2018
    --> UPDATED_GFED4          :       true     # 2009-2019
    --> JPL_WETCHARTS          :       true     # 2009-2017
    --> SEEPS                  :       true     # 2012
    --> LAKES                  :       false    # 2009-2015
    --> FUNG_TERMITES          :       true     # 1985
    --> FUNG_SOIL_ABSORPTION   :       false    # 2009-2015
    --> MeMo_SOIL_ABSORPTION   :       true     # 1990-2009 or clim.
    --> CMIP6_SFC_LAND_ANTHRO  :       false    # 1850-2100
    --> CMIP6_SHIP             :       false    # 1850-2100
    --> BB4MIPS                :       false    # 1850-2100
# ..... Non-Emissions Data ...........
    --> GMI_CH4_LOSS           :       true     # 1985
    --> GLOBAL_CL              :       true     # 2010-2019
# ..... Options for the IMI ..........
    --> AnalyticalInv          :       false
    --> Emis_ScaleFactor       :       false
    --> OH_ScaleFactor         :       false
# ----- CO and CO2-only INVENTORIES AND DATA ----------------------------------
# ..... Global inventories ...........
    --> AEIC2019_DAILY         :       false    # 2019 (daily data)
    --> AEIC2019_MONMEAN       :       true     # 2019 (monthly-mean data)
    --> AEIC_SCALE_1990_2019   :       true     # Scale to year in 1990-2019
    --> CEDSv2                 :       true     # 1750-2019
    --> CEDS_GBDMAPS           :       false    # 1970-2017
    --> CEDS_GBDMAPS_byFuelType:       false    # 1970-2017
    --> HTAP                   :       false    # 2008-2010
    --> FOSSIL_ODIAC           :       true     # 2000-2018
    --> FOSSIL_CDIAC           :       false    # 1980-2014
    --> OCEAN_EXCH_TAKA09      :       false    # 2000
    --> OCEAN_EXCH_SCALED      :       true     # 2000-2013
    --> BBIO_DIURNAL           :       false    # 1985
    --> BBIO_SIB3              :       true     # 2006-2010
    --> NET_TERR_EXCH          :       true     # 2000
    --> CO2CORR                :       true     # 2000-2018
# ..... Regional inventories .........
    --> APEI                   :       false    # 1989-2014
    --> NEI2016_MONMEAN        :       false    # 2002-2020
    --> DICE_Africa            :       false    # 2013
# ..... Non-emissions data ...........
    --> CO2_COPROD             :       true     # 2004-2009
    --> PROD_CO_CH4            :       true     # 2010-2019
    --> PROD_CO_NMVOC          :       true     # 2010-2019
    --> GMI_PROD_LOSS          :       true     # 2005
    --> GMD_SFC_CH4            :       true     # 1979-2020
    --> CMIP6_SFC_CH4          :       false    # 1750-1978
# ..... Ship emissions ...............
    --> CEDSv2_SHIP            :       true     # 1750-2017
    --> CEDS_GBDMAPS_SHIP      :       false    # 1970-2017
    --> CEDS_GBDMAPS_SHIP_byFuelType:  false    # 1970-2017
    --> ICOADS_SHIP            :       false    # 2004
#------ OCS FLUX DATA ---------------------------------------------------------
    --> OCS_ANTHRO_FLUX        :       true     # 2012-2016
    --> OCS_BIOMASS_FLUX       :       true     # 2012-2016
    --> OCS_MISSING_OCEAN_FLUX :       true     # 2012-2016
    --> OCS_OCEAN_FLUX         :       true     # 2012-2016
# ----- FUTURE EMISSIONS ------------------------------------------------------
    --> RCP_3PD                :       false    # 2005-2100
    --> RCP_45                 :       false    # 2005-2100
    --> RCP_60                 :       false    # 2005-2100
    --> RCP_85                 :       false    # 2005-2100
# ----- NON-EMISSIONS DATA ----------------------------------------------------
    --> OLSON_LANDMAP          :       true     # 1985
    --> YUAN_MODIS_LAI         :       true     # 2000-2020
    --> GLOBAL_OH              :       true     # 2010-2019
    --> GLOBAL_OH_GCv5         :       false    # 1985 (GCv5 OH, need for IMI)
# -----------------------------------------------------------------------------
111     GFED                   : on    CH4/CO/CO2
    --> GFED4                  :       true
    --> GFED_daily             :       true
    --> GFED_3hourly           :       false
    --> Scaling_CO             :       1.05

### END SECTION EXTENSION SWITCHES ###

###############################################################################
### BEGIN SECTION BASE EMISSIONS
###############################################################################

# ExtNr	Name sourceFile	sourceVar sourceTime C/R/E SrcDim SrcUnit Species ScalIDs Cat Hier

(((EMISSIONS

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%                        CH4 EMISSIONS                               %%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
(((USE_CH4_DATA

#==============================================================================
# --- CH4: Gridded EPA (Maasakkers et al., Environ. Sci. Technol., 2016) ---
#
# NOTES:
# - Use Hier=100 to add to Canada and Mexico regional inventories
# - Do not use GEPA forest fire emissions. Use QFED instead.
# - Apply seasonal scale factors to manure and rice emissions
#==============================================================================
(((GEPA
0 GEPA_OIL              $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_1B2a_Petroleum                     2012/1/1/0    C xy molec/cm2/s CH4 1008    1 100
0 GEPA_GAS_PRODUCTION   $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_1B2b_Natural_Gas_Production        2012/1/1/0    C xy molec/cm2/s CH4 1008    2 100
0 GEPA_GAS_PROCESSING   $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_1B2b_Natural_Gas_Processing        2012/1/1/0    C xy molec/cm2/s CH4 1008    2 100
0 GEPA_GAS_TRANSMISSION $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_1B2b_Natural_Gas_Transmission      2012/1/1/0    C xy molec/cm2/s CH4 1008    2 100
0 GEPA_GAS_DISTRIBUTION $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_1B2b_Natural_Gas_Distribution      2012/1/1/0    C xy molec/cm2/s CH4 1008    2 100
0 GEPA_COAL_UNDERGROUND $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_1B1a_Coal_Mining_Underground       2012/1/1/0    C xy molec/cm2/s CH4 1008    3 100
0 GEPA_COAL_SURFACE     $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_1B1a_Coal_Mining_Surface           2012/1/1/0    C xy molec/cm2/s CH4 1008    3 100
0 GEPA_COAL_ABANDONED   $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_1B1a_Abandoned_Coal                2012/1/1/0    C xy molec/cm2/s CH4 1008    3 100
0 GEPA_LIVESTOCK__4A    $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_4A_Enteric_Fermentation            2012/1/1/0    C xy molec/cm2/s CH4 1008    4 100
0 GEPA_LIVESTOCK__4B    $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_4B_Manure_Management               2012/1/1/0    C xy molec/cm2/s CH4 10/1008 4 100
0 GEPA_LANDFILLS_MUNI   $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_6A_Landfills_Municipal             2012/1/1/0    C xy molec/cm2/s CH4 1008    5 100
0 GEPA_LANDFILLS_INDU   $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_6A_Landfills_Industrial            2012/1/1/0    C xy molec/cm2/s CH4 1008    5 100
0 GEPA_WASTEWATER_DOME  $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_6B_Wastewater_Treatment_Domestic   2012/1/1/0    C xy molec/cm2/s CH4 1008    6 100
0 GEPA_WASTEWATER_INDU  $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_6B_Wastewater_Treatment_Industrial 2012/1/1/0    C xy molec/cm2/s CH4 1008    6 100
0 GEPA_RICE             $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc emissions_4C_Rice_Cultivation                 2012/1/1/0    C xy molec/cm2/s CH4 11/1008 7 100
0 GEPA_OTHER__1A_M      $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_1A_Combustion_Mobile               2012/1/1/0    C xy molec/cm2/s CH4 1008    8 100
0 GEPA_OTHER__1A_S      $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_1A_Combustion_Stationary           2012/1/1/0    C xy molec/cm2/s CH4 1008    8 100
0 GEPA_OTHER__2B5       $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_2B5_Petrochemical_Production       2012/1/1/0    C xy molec/cm2/s CH4 1008    8 100
0 GEPA_OTHER__2C2       $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_2C2_Ferroalloy_Production          2012/1/1/0    C xy molec/cm2/s CH4 1008    8 100
0 GEPA_OTHER__4F        $ROOT/CH4/v2022-11/GEPA/GEPA_Monthly.nc emissions_4F_Field_Burning                   2012/1-12/1/0 C xy molec/cm2/s CH4 1008    8 100
#0 GEPA_OTHER__5        $ROOT/CH4/v2022-11/GEPA_Daily.nc        emissions_5_Forest_Fires                     2012/1/1/0    C xy molec/cm2/s CH4 1008    9 100
0 GEPA_OTHER__6D        $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_6D_Composting                      2012/1/1/0    C xy molec/cm2/s CH4 1008    8 100
0 GEPA_COAST_OIL        $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_1B2a_Petroleum                     2012/1/1/0    C xy molec/cm2/s CH4 1009    1 5
0 GEPA_COAST_GAS_PR     $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_1B2b_Natural_Gas_Production        2012/1/1/0    C xy molec/cm2/s CH4 1009    2 5
0 GEPA_COAST_GAS_PC     $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_1B2b_Natural_Gas_Processing        2012/1/1/0    C xy molec/cm2/s CH4 1009    2 5
0 GEPA_COAST_GAS_TR     $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_1B2b_Natural_Gas_Transmission      2012/1/1/0    C xy molec/cm2/s CH4 1009    2 5
0 GEPA_COAST_GAS_DS     $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_1B2b_Natural_Gas_Distribution      2012/1/1/0    C xy molec/cm2/s CH4 1009    2 5
0 GEPA_COAST_COAL_U     $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_1B1a_Coal_Mining_Underground       2012/1/1/0    C xy molec/cm2/s CH4 1009    3 5
0 GEPA_COAST_COAL_S     $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_1B1a_Coal_Mining_Surface           2012/1/1/0    C xy molec/cm2/s CH4 1009    3 5
0 GEPA_COAST_COAL_A     $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_1B1a_Abandoned_Coal                2012/1/1/0    C xy molec/cm2/s CH4 1009    3 5
0 GEPA_COAST_LVSTK_F    $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_4A_Enteric_Fermentation            2012/1/1/0    C xy molec/cm2/s CH4 1009    4 5
0 GEPA_COAST_LVSTK_M    $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_4B_Manure_Management               2012/1/1/0    C xy molec/cm2/s CH4 10/1009 4 5
0 GEPA_COAST_LDF_M      $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_6A_Landfills_Municipal             2012/1/1/0    C xy molec/cm2/s CH4 1009    5 5
0 GEPA_COAST_LDF_I      $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_6A_Landfills_Industrial            2012/1/1/0    C xy molec/cm2/s CH4 1009    5 5
0 GEPA_COAST_WSTW_D     $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_6B_Wastewater_Treatment_Domestic   2012/1/1/0    C xy molec/cm2/s CH4 1009    6 5
0 GEPA_COAST_WSTW_I     $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_6B_Wastewater_Treatment_Industrial 2012/1/1/0    C xy molec/cm2/s CH4 1009    6 5
0 GEPA_COAST_RICE       $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_4C_Rice_Cultivation                2012/1/1/0    C xy molec/cm2/s CH4 11/1009 7 5
0 GEPA_COAST_OTH_1AM    $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_1A_Combustion_Mobile               2012/1/1/0    C xy molec/cm2/s CH4 1009    8 5
0 GEPA_COAST_OTH_1AS    $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_1A_Combustion_Stationary           2012/1/1/0    C xy molec/cm2/s CH4 1009    8 5
0 GEPA_COAST_OTH_2B5    $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_2B5_Petrochemical_Production       2012/1/1/0    C xy molec/cm2/s CH4 1009    8 5
0 GEPA_COAST_OTH_2C2    $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_2C2_Ferroalloy_Production          2012/1/1/0    C xy molec/cm2/s CH4 1009    8 5
0 GEPA_COAST_OTH_4F     $ROOT/CH4/v2022-11/GEPA/GEPA_Monthly.nc emissions_4F_Field_Burning                   2012/1-12/1/0 C xy molec/cm2/s CH4 1009    8 5
#0 GEPA_COAST_OTH_5     $ROOT/CH4/v2022-11/GEPA_Daily.nc        emissions_5_Forest_Fires                     2012/1/1/0    C xy molec/cm2/s CH4 1009    9 5
0 GEPA_COAST_OTH_6D     $ROOT/CH4/v2022-11/GEPA/GEPA_Annual.nc  emissions_6D_Composting                      2012/1/1/0    C xy molec/cm2/s CH4 1009    8 5
)))GEPA

#==============================================================================
# --- CH4: Mexico emissions (Scarpelli et. al, Environ. Res. Lett., 2020) ---
#
# NOTES:
# - Use Hier=100 to add to Canada and USA regional inventories
#==============================================================================
(((Scarpelli_Mexico
0 MEX_OIL               $ROOT/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_oil_2015.nc          emis_ch4 2015/1/1/0 C xy molec/cm2/s CH4 1001    1 100
0 MEX_GAS               $ROOT/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_gas_2015.nc          emis_ch4 2015/1/1/0 C xy molec/cm2/s CH4 1001    2 100
0 MEX_COAL              $ROOT/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_coal_2015.nc         emis_ch4 2015/1/1/0 C xy molec/cm2/s CH4 1001    3 100
0 MEX_LIVESTOCK_A       $ROOT/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_livestock_A_2015.nc  emis_ch4 2015/1/1/0 C xy molec/cm2/s CH4 1001    4 100
0 MEX_LIVESTOCK_B       $ROOT/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_livestock_B_2015.nc  emis_ch4 2015/1/1/0 C xy molec/cm2/s CH4 10/1001 4 100
0 MEX_LANDFILLS         $ROOT/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_landfill_2015.nc     emis_ch4 2015/1/1/0 C xy molec/cm2/s CH4 1001    5 100
0 MEX_WASTEWATER        $ROOT/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_waste_2015.nc        emis_ch4 2015/1/1/0 C xy molec/cm2/s CH4 1001    6 100
0 MEX_RICE              $ROOT/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_rice_2015.nc         emis_ch4 2015/1/1/0 C xy molec/cm2/s CH4 11/1001 7 100
0 MEX_OTHER             $ROOT/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_other_anthro_2015.nc emis_ch4 2015/1/1/0 C xy molec/cm2/s CH4 1001    8 100

### Make sure to include offshore/coastal emissions (Hier=1 to add to EDGAR, Hier=5 to add to GFEI) ###
0 MEX_OIL_COAST         $ROOT/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_oil_2015.nc          emis_ch4 2015/1/1/0 C xy molec/cm2/s CH4 1010    1 5
0 MEX_GAS_COAST         $ROOT/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_gas_2015.nc          emis_ch4 2015/1/1/0 C xy molec/cm2/s CH4 1010    2 5
0 MEX_COAL_COAST        $ROOT/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_coal_2015.nc         emis_ch4 2015/1/1/0 C xy molec/cm2/s CH4 1010    3 5
0 MEX_LIVESTOCK_A_COAST $ROOT/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_livestock_A_2015.nc  emis_ch4 2015/1/1/0 C xy molec/cm2/s CH4 1010    4 1
0 MEX_LIVESTOCK_B_COAST $ROOT/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_livestock_B_2015.nc  emis_ch4 2015/1/1/0 C xy molec/cm2/s CH4 10/1010 4 1
0 MEX_LANDFILLS_COAST   $ROOT/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_landfill_2015.nc     emis_ch4 2015/1/1/0 C xy molec/cm2/s CH4 1010    5 1
0 MEX_WASTEWATER_COAST  $ROOT/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_waste_2015.nc        emis_ch4 2015/1/1/0 C xy molec/cm2/s CH4 1010    6 1
0 MEX_RICE_COAST        $ROOT/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_rice_2015.nc         emis_ch4 2015/1/1/0 C xy molec/cm2/s CH4 11/1010 7 1
0 MEX_OTHER_COAST       $ROOT/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_other_anthro_2015.nc emis_ch4 2015/1/1/0 C xy molec/cm2/s CH4 1010    8 1
)))Scarpelli_Mexico

#==============================================================================
# --- CH4: Canada emissions (Scarpelli et al., Environ. Res. Lett., 2022) ---
#
# NOTES:
# - Use Hier=100 to add to USA and Mexico regional inventories
#==============================================================================
(((Scarpelli_Canada
0 CAN_OIL_GAS_COMBUSTION  $ROOT/CH4/v2022-01/Scarpelli_Canada/can_emis_oil_gas_combustion_2018.nc  oil_gas_combustion_total  2018/1/1/0 C xy kg/m2/s CH4 1002 1/2 100
0 CAN_OIL_GAS_LEAKAGE     $ROOT/CH4/v2022-01/Scarpelli_Canada/can_emis_oil_gas_leakage_2018.nc     oil_gas_leakage_total     2018/1/1/0 C xy kg/m2/s CH4 1002 1/2 100
0 CAN_OIL_GAS_VENT_FLARE  $ROOT/CH4/v2022-01/Scarpelli_Canada/can_emis_oil_gas_vent_flare_2018.nc  oil_gas_vent_flare_total  2018/1/1/0 C xy kg/m2/s CH4 1002 1/2 100
0 CAN_COAL                $ROOT/CH4/v2022-01/Scarpelli_Canada/can_emis_coal_2018.nc                coal_total                2018/1/1/0 C xy kg/m2/s CH4 1002 3   100
0 CAN_LIVESTOCK           $ROOT/CH4/v2022-01/Scarpelli_Canada/can_emis_livestock_2018.nc           livestock_total           2018/1/1/0 C xy kg/m2/s CH4 1002 4   100
0 CAN_SOLID_WASTE         $ROOT/CH4/v2022-01/Scarpelli_Canada/can_emis_solid_waste_2018.nc         solid_waste_total         2018/1/1/0 C xy kg/m2/s CH4 1002 5   100
0 CAN_WASTEWATER          $ROOT/CH4/v2022-01/Scarpelli_Canada/can_emis_wastewater_2018.nc          wastewater_total          2018/1/1/0 C xy kg/m2/s CH4 1002 6   100
0 CAN_OTHER               $ROOT/CH4/v2022-01/Scarpelli_Canada/can_emis_other_minor_sources_2018.nc other_minor_sources_total 2018/1/1/0 C xy kg/m2/s CH4 1002 8   100

### Make sure to include offshore/coastal emissions (Hier=1 to add to EDGAR, Hier=5 to add to GFEI) ###
0 CAN_OIL_GAS_COMBUSTION_COAST  $ROOT/CH4/v2022-01/Scarpelli_Canada/can_emis_oil_gas_combustion_2018.nc  oil_gas_combustion_total  2018/1/1/0 C xy kg/m2/s CH4 1011 1/2 5
0 CAN_OIL_GAS_LEAKAGE_COAST     $ROOT/CH4/v2022-01/Scarpelli_Canada/can_emis_oil_gas_leakage_2018.nc     oil_gas_leakage_total     2018/1/1/0 C xy kg/m2/s CH4 1011 1/2 5
0 CAN_OIL_GAS_VENT_FLARE_COAST  $ROOT/CH4/v2022-01/Scarpelli_Canada/can_emis_oil_gas_vent_flare_2018.nc  oil_gas_vent_flare_total  2018/1/1/0 C xy kg/m2/s CH4 1011 1/2 5
0 CAN_COAL_COAST                $ROOT/CH4/v2022-01/Scarpelli_Canada/can_emis_coal_2018.nc                coal_total                2018/1/1/0 C xy kg/m2/s CH4 1011 3   5
0 CAN_LIVESTOCK_COAST           $ROOT/CH4/v2022-01/Scarpelli_Canada/can_emis_livestock_2018.nc           livestock_total           2018/1/1/0 C xy kg/m2/s CH4 1011 4   1
0 CAN_SOLID_WASTE_COAST         $ROOT/CH4/v2022-01/Scarpelli_Canada/can_emis_solid_waste_2018.nc         solid_waste_total         2018/1/1/0 C xy kg/m2/s CH4 1011 5   1
0 CAN_WASTEWATER_COAST          $ROOT/CH4/v2022-01/Scarpelli_Canada/can_emis_wastewater_2018.nc          wastewater_total          2018/1/1/0 C xy kg/m2/s CH4 1011 6   1
0 CAN_OTHER_COAST               $ROOT/CH4/v2022-01/Scarpelli_Canada/can_emis_other_minor_sources_2018.nc other_minor_sources_total 2018/1/1/0 C xy kg/m2/s CH4 1011 8   1
)))Scarpelli_Canada

#==============================================================================
# --- CH4: Global Fuel Exploitation Inventory (GFEI v2, Scarpelli et al., 2021)
#
# This inventory will replace EDGAR (oil, gas, & coal)
#==============================================================================
(((GFEIv2
0 GFEI_CH4_OIL  $ROOT/CH4/v2022-01/GFEIv2/Global_Fuel_Exploitation_Inventory_v2_2019_Oil_All.nc  emis_ch4 2019/1/1/0 C xy molec/cm2/s CH4 -  1 5
0 GFEI_CH4_GAS  $ROOT/CH4/v2022-01/GFEIv2/Global_Fuel_Exploitation_Inventory_v2_2019_Gas_All.nc  emis_ch4 2019/1/1/0 C xy molec/cm2/s CH4 -  2 5
0 GFEI_CH4_COAL $ROOT/CH4/v2022-01/GFEIv2/Global_Fuel_Exploitation_Inventory_v2_2019_Coal.nc     emis_ch4 2019/1/1/0 C xy molec/cm2/s CH4 -  3 5
)))GFEIv2

#==============================================================================
# --- CH4: EDGAR v6.0 emissions, various sectors ---
#==============================================================================
(((EDGARv6
0 CH4_OIL__1B2a          $ROOT/CH4/v2022-11/EDGARv6/$YYYY/v6.0_CH4_$YYYY_PRO_OIL.0.1x0.1.nc          emi_ch4 2000-2018/1-12/1/0 C xy kg/m2/s CH4 - 1 1
0 CH4_OIL__1B2c          $ROOT/CH4/v2022-11/EDGARv6/$YYYY/v6.0_CH4_$YYYY_PRO_GAS.0.1x0.1.nc          emi_ch4 2000-2018/1-12/1/0 C xy kg/m2/s CH4 - 2 1
0 CH4_COAL__1B1a         $ROOT/CH4/v2022-11/EDGARv6/$YYYY/v6.0_CH4_$YYYY_PRO_COAL.0.1x0.1.nc         emi_ch4 2000-2018/1-12/1/0 C xy kg/m2/s CH4 - 3 1
0 CH4_LIVESTOCK__4A      $ROOT/CH4/v2022-11/EDGARv6/$YYYY/v6.0_CH4_$YYYY_ENF.0.1x0.1.nc              emi_ch4 2000-2018/1-12/1/0 C xy kg/m2/s CH4 - 4 1
0 CH4_LIVESTOCK__4B      $ROOT/CH4/v2022-11/EDGARv6/$YYYY/v6.0_CH4_$YYYY_MNM.0.1x0.1.nc              emi_ch4 2000-2018/1-12/1/0 C xy kg/m2/s CH4 - 4 1
0 CH4_LANDFILLS__6A_6D   $ROOT/CH4/v2022-11/EDGARv6/$YYYY/v6.0_CH4_$YYYY_SWD_LDF.0.1x0.1.nc          emi_ch4 2000-2018/1-12/1/0 C xy kg/m2/s CH4 - 5 1
0 CH4_WASTEWATER__6B     $ROOT/CH4/v2022-11/EDGARv6/$YYYY/v6.0_CH4_$YYYY_WWT.0.1x0.1.nc              emi_ch4 2000-2018/1-12/1/0 C xy kg/m2/s CH4 - 6 1
0 CH4_RICE__4C_4D        $ROOT/CH4/v2022-11/EDGARv6/$YYYY/v6.0_CH4_$YYYY_AGS.0.1x0.1.nc              emi_ch4 2000-2018/1-12/1/0 C xy kg/m2/s CH4 - 7 1
0 CH4_OTHER__1A1_1B1_1B2 $ROOT/CH4/v2022-11/EDGARv6/$YYYY/v6.0_CH4_$YYYY_REF_TRF.0.1x0.1.nc          emi_ch4 2000-2018/1-12/1/0 C xy kg/m2/s CH4 - 8 1
0 CH4_OTHER__1A1a        $ROOT/CH4/v2022-11/EDGARv6/$YYYY/v6.0_CH4_$YYYY_ENE.0.1x0.1.nc              emi_ch4 2000-2018/1-12/1/0 C xy kg/m2/s CH4 - 8 1
0 CH4_OTHER__1A2         $ROOT/CH4/v2022-11/EDGARv6/$YYYY/v6.0_CH4_$YYYY_IND.0.1x0.1.nc              emi_ch4 2000-2018/1-12/1/0 C xy kg/m2/s CH4 - 8 1
0 CH4_OTHER__1A3a_CDS    $ROOT/CH4/v2022-11/EDGARv6/$YYYY/v6.0_CH4_$YYYY_TNR_Aviation_CDS.0.1x0.1.nc emi_ch4 2000-2018/1-12/1/0 C xy kg/m2/s CH4 - 8 1
0 CH4_OTHER__1A3a_CRS    $ROOT/CH4/v2022-11/EDGARv6/$YYYY/v6.0_CH4_$YYYY_TNR_Aviation_CRS.0.1x0.1.nc emi_ch4 2000-2018/1-12/1/0 C xy kg/m2/s CH4 - 8 1
0 CH4_OTHER__1A3a_LTO    $ROOT/CH4/v2022-11/EDGARv6/$YYYY/v6.0_CH4_$YYYY_TNR_Aviation_LTO.0.1x0.1.nc emi_ch4 2000-2018/1-12/1/0 C xy kg/m2/s CH4 - 8 1
0 CH4_OTHER__1A3b        $ROOT/CH4/v2022-11/EDGARv6/$YYYY/v6.0_CH4_$YYYY_TRO_noRES.0.1x0.1.nc        emi_ch4 2000-2018/1-12/1/0 C xy kg/m2/s CH4 - 8 1
0 CH4_OTHER__1A3c_1A3e   $ROOT/CH4/v2022-11/EDGARv6/$YYYY/v6.0_CH4_$YYYY_TNR_Other.0.1x0.1.nc        emi_ch4 2000-2018/1-12/1/0 C xy kg/m2/s CH4 - 8 1
0 CH4_OTHER__1A3d_1C2    $ROOT/CH4/v2022-11/EDGARv6/$YYYY/v6.0_CH4_$YYYY_TNR_Ship.0.1x0.1.nc         emi_ch4 2000-2018/1-12/1/0 C xy kg/m2/s CH4 - 8 1
0 CH4_OTHER__1A4         $ROOT/CH4/v2022-11/EDGARv6/$YYYY/v6.0_CH4_$YYYY_RCO.0.1x0.1.nc              emi_ch4 2000-2018/1-12/1/0 C xy kg/m2/s CH4 - 8 1
0 CH4_OTHER__2B          $ROOT/CH4/v2022-11/EDGARv6/$YYYY/v6.0_CH4_$YYYY_CHE.0.1x0.1.nc              emi_ch4 2000-2018/1-12/1/0 C xy kg/m2/s CH4 - 8 1
0 CH4_OTHER__2C          $ROOT/CH4/v2022-11/EDGARv6/$YYYY/v6.0_CH4_$YYYY_IRO.0.1x0.1.nc              emi_ch4 2000-2018/1-12/1/0 C xy kg/m2/s CH4 - 8 1
0 CH4_OTHER__4F          $ROOT/CH4/v2022-11/EDGARv6/$YYYY/v6.0_CH4_$YYYY_AWB.0.1x0.1.nc              emi_ch4 2000-2018/1-12/1/0 C xy kg/m2/s CH4 - 8 1
0 CH4_OTHER__6C          $ROOT/CH4/v2022-11/EDGARv6/$YYYY/v6.0_CH4_$YYYY_SWD_INC.0.1x0.1.nc          emi_ch4 2000-2018/1-12/1/0 C xy kg/m2/s CH4 - 8 1
)))EDGARv6

#==============================================================================
# --- CH4: CEDS (historical) or Shared Socioeconomic Pathways (future) ---
#==============================================================================
(((CMIP6_SFC_LAND_ANTHRO
0 CMIP6_CH4_AGR   $ROOT/CMIP6/v2021-01/$GCAP2SCENARIO/$GCAP2SCENARIO_$YYYY.nc4 CH4_agr 1970-2014/1-12/1/0 C xy kg/m2/s CH4 - 4 1
0 CMIP6_CH4_ENE   $ROOT/CMIP6/v2021-01/$GCAP2SCENARIO/$GCAP2SCENARIO_$YYYY.nc4 CH4_ene 1970-2014/1-12/1/0 C xy kg/m2/s CH4 - 2 1
0 CMIP6_CH4_IND   $ROOT/CMIP6/v2021-01/$GCAP2SCENARIO/$GCAP2SCENARIO_$YYYY.nc4 CH4_ind 1970-2014/1-12/1/0 C xy kg/m2/s CH4 - 8 1
0 CMIP6_CH4_TRA   $ROOT/CMIP6/v2021-01/$GCAP2SCENARIO/$GCAP2SCENARIO_$YYYY.nc4 CH4_tra 1970-2014/1-12/1/0 C xy kg/m2/s CH4 - 8 1
0 CMIP6_CH4_RCO   $ROOT/CMIP6/v2021-01/$GCAP2SCENARIO/$GCAP2SCENARIO_$YYYY.nc4 CH4_rco 1970-2014/1-12/1/0 C xy kg/m2/s CH4 - 2 1
0 CMIP6_CH4_SLV   $ROOT/CMIP6/v2021-01/$GCAP2SCENARIO/$GCAP2SCENARIO_$YYYY.nc4 CH4_slv 1970-2014/1-12/1/0 C xy kg/m2/s CH4 - 8 1
0 CMIP6_CH4_WST   $ROOT/CMIP6/v2021-01/$GCAP2SCENARIO/$GCAP2SCENARIO_$YYYY.nc4 CH4_wst 1970-2014/1-12/1/0 C xy kg/m2/s CH4 - 5 1
)))CMIP6_SFC_LAND_ANTHRO
(((CMIP6_SHIP
0 CMIP6_CH4_SHP   $ROOT/CMIP6/v2021-01/$GCAP2SCENARIO/$GCAP2SCENARIO_$YYYY.nc4 CH4_shp 1970-2014/1-12/1/0 C xy kg/m2/s CH4 - 8 1
)))CMIP6_SHIP

#==============================================================================
# --- CH4: BB4MIPs historical / SSP future biomass burning inventories ---
#==============================================================================
(((BB4MIPS
0 CMIP6_BB_CH4    $ROOT/CMIP6/v2021-01/$GCAP2SCENARIO/$GCAP2SCENARIO_$YYYY.nc4 CH4_bbn 1750-2015/1-12/1/0 C xyL=1:PBL kg/m2/s CH4 75 9 1
)))BB4MIPS

#==============================================================================
# --- CH4: biomass burning ---
#==============================================================================

# --- QFED v2.4r8 ---
# NOTE: Use this instead of data in QFED/v2014-09 (these are buggy)
(((QFED2
0 QFED_CH4   $ROOT/CH4/v2017-10/QFED/QFEDv2.4r8.emis_ch4.$YYYY.nc              emi_ch4 2009-2015/1-12/1-31/0 C xy molec/cm2/s CH4 -  9 2
)))QFED2

# --- Updated GFED4 biomass burning ---
# Emissions thru 2019 (cf Yuzhong Zhang)
(((UPDATED_GFED4
0 UPDATED_GFED4_CH4    $ROOT/CH4/v2020-09/GFED/GFED4_CH4_FIRE_$YYYY.nc emi_ch4 2009-2019/1-12/1-31/0 C xy molec/cm2/s CH4 -  9 3
)))UPDATED_GFED4

#==============================================================================
# --- CH4: JPL WetCHARTs v1.0 (Bloom et al., https://doi.org/10.3334/ORNLDAAC/1502) ---
#==============================================================================
(((JPL_WETCHARTS
0 JPLW_CH4  $ROOT/CH4/v2020-09/JPL_WetCharts/HEensemble/JPL_WetCharts_2010-2019.Ensemble_Mean.0.5x0.5.nc emi_ch4 2010-2019/1-12/1/0 C xy molec/cm2/s CH4 - 10 1
)))JPL_WETCHARTS

#==============================================================================
# --- CH4: Geological Seeps ---
#==============================================================================
(((SEEPS
0 CH4_SEEPS $ROOT/CH4/v2020-04/Seeps/Etiope_CH4GeologicalEmis_ScaledToHmiel.1x1.nc emi_ch4 2012/1/1/0 C xy kg/m2/s CH4 - 11 1
)))SEEPS

#==============================================================================
# --- CH4: Emissions from Lakes (Maasakkers et al., 2019) ---
#==============================================================================
(((LAKES
0 CH4_LAKES $ROOT/CH4/v2017-10/Lakes/Maasakkers_Lakes_$YYYY.01x01.nc emi_ch4 $YYYY/1-12/1/0 C xy molec/cm2/s CH4 - 12 1
)))LAKES

#==============================================================================
# --- CH4: Soil absorption & termites from Fung et al, 1991 ---
#
# NOTES:
# - Multiply soil absorption by -1 to get a "negative" flux.
#   (Only apply the scaling factor when adding to the total CH4 simulation)
# - Use updated soil absorption emissions with seasonality applied
# - For tagged CH4 simulations soil absorption is treated as a loss and applied
#   in mixing_mod.F
#==============================================================================
(((FUNG_TERMITES
0 CH4_TERMITES     $ROOT/CH4/v2022-11/4x5/termites.geos.4x5.nc CH4 1985/1/1/0 C xy kg/m2/s CH4 - 13 1
)))FUNG_TERMITES
(((FUNG_SOIL_ABSORPTION
0 CH4_SOILABSORB   $ROOT/CH4/v2019-10/Fung_SoilAbs/Soil_Absorption_4x5_$YYYY.nc CH4 2009-2015/1-12/1/0 C xy molec/cm2/s CH4 1 14 1
)))FUNG_SOIL_ABSORPTION

#==============================================================================
# --- CH4: Soil absorption from MeMo model (Murguia-Flores et al. 2018, GMD) ---
#
# - Multiply soil absorption by -1 to get a "negative" flux
#==============================================================================
(((MeMo_SOIL_ABSORPTION
# Climatology
0 CH4_SOILABSORB $ROOT/CH4/v2019-10/MeMo_SoilAbs/MeMo_CH4uptake_Climatology.nc CH4uptake 2009/1-12/1/0 C xy kg/m2/s CH4 1 14 2

# 1990-2009
#0 CH4_SOILABSORB $ROOT/CH4/v2019-10/MeMo_SoilAbs/MeMo_CH4uptake_1990-2009.nc CH4uptake 1990-2009/1-12/1/0 C xy kg/m2/s CH4 1 14 2
)))MeMo_SOIL_ABSORPTION

#==============================================================================
# --- CH4: RCP future emission scenarios ---
#==============================================================================
(((RCP_3PD
0 RCP3PD_CH4 $ROOT/RCP/v2020-07/RCP_3PD/RCPs_anthro_CH4_2005-2100_23474.nc  ACCMIP 2005-2100/1/1/0 ID xy kg/m2/s  CH4 - 1 1
)))RCP_3PD

(((RCP_45
0 RCP45_CH4  $ROOT/RCP/v2020-07/RCP_45/RCPs_anthro_CH4_2005-2100_27424.nc   ACCMIP 2005-2100/1/1/0 ID xy kg/m2/s  CH4 - 1 1
)))RCP_45

(((RCP_60
0 RCP60_CH4  $ROOT/RCP/v2020-07/RCP_60/RCPs_anthro_CH4_2005-2100_43190.nc   ACCMIP 2005-2100/1/1/0 ID xy kg/m2/s  CH4 - 1 1
)))RCP_60

(((RCP_85
0 RCP85_CH4  $ROOT/RCP/v2020-07/RCP_85/RCPs_anthro_CH4_2005-2100_43533.nc   ACCMIP 2005-2100/1/1/0 ID xy kg/m2/s  CH4 - 1 1
)))RCP_85

)))USE_CH4_DATA

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%                        CO EMISSIONS                                %%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
(((USE_CO_DATA

#==============================================================================
# --- CO: APEI (Canada) ---
#==============================================================================
(((APEI
0 APEI_CO    $ROOT/APEI/v2016-11/APEI.0.1x0.1.nc CO 1989-2014/1/1/0 RF xy kg/m2/s CO    26/52/2002      1 30
)))APEI

#==============================================================================
# --- CO: EPA NEI2016 v1 (USA) ---
#
# NOTES:
#  * Barron Henderson wrote, "The EPA emission modeling platform always
#    includes our best estimate of that year's emissions for Canada and Mexico
#    (othpt, othar, ptfire_oth). Using that estimate is likely good. However,
#    that can lead to a discontinuity in Mexico and Canada."
#    - By default only emissions over the CONUS are used (via Mask #2007)
#    - To include emissions over Canada and Mexico, users may revert to the old
#      US mask file $ROOT/MASKS/v2018-09/USA_LANDMASK_NEI2011_0.1x0.1.20160921.nc
#  * The base year of these emissions is 2016 and emissions are scaled to
#    2002-2020 using data from the EPA Trends Report for Tier 1 CAPS (obtained
#    21 Sep 2021).
#    - See NEI2016/v2021-06/national_tier1_caps+HEMCOscaling.xlsx for details.
#    - Cl2 and HCl emissions are only available for 2016 because those species
#      aren't included in the Tier 1 CAPS file/
#==============================================================================
(((NEI2016_MONMEAN
0 EPA16_CO__airportsCO   $ROOT/NEI2016/v2021-06/2016fh_16j_airports_0pt1degree_month_$MM.ncf                   CO  2002-2020/1-12/1/0 RF xyz  kg/m2/s CO 26/211/252/2007 1 50
0 EPA16_CO__nonptCO      $ROOT/NEI2016/v2021-06/2016fh_16j_nonpt_0pt1degree_month_$MM.ncf                      CO  2002-2020/1-12/1/0 RF xyz  kg/m2/s CO 26/211/252/2007 1 50
0 EPA16_CO__nonroadCO    $ROOT/NEI2016/v2021-06/2016fh_16j_nonroad_0pt1degree_month_$MM.ncf                    CO  2002-2020/1-12/1/0 RF xyz  kg/m2/s CO 26/211/252/2007 1 50
0 EPA16_CO__npogCO       $ROOT/NEI2016/v2021-06/2016fh_16j_np_oilgas_0pt1degree_month_$MM.ncf                  CO  2002-2020/1-12/1/0 RF xyz  kg/m2/s CO 26/211/252/2007 1 50
0 EPA16_CO__onroadCO     $ROOT/NEI2016/v2021-06/2016fh_16j_onroad_0pt1degree_month_$MM.ncf                     CO  2002-2020/1-12/1/0 RF xyz  kg/m2/s CO 26/211/252/2007 1 50
0 EPA16_CO__onroad_caCO  $ROOT/NEI2016/v2021-06/2016fh_16j_onroad_ca_adj_0pt1degree_month_$MM.ncf              CO  2002-2020/1-12/1/0 RF xyz  kg/m2/s CO 26/211/252/2007 1 50
0 EPA16_CO__railCO       $ROOT/NEI2016/v2021-06/2016fh_16j_rail_0pt1degree_month_$MM.ncf                       CO  2002-2020/1-12/1/0 RF xyz  kg/m2/s CO 26/211/252/2007 1 50
0 EPA16_CO__c1c2CO       $ROOT/NEI2016/v2021-06/2016fh_16j_emln_cmv_c1c2_12_0pt1degree_3D_month_$MM.ncf        CO  2002-2020/1-12/1/0 RF xyz  kg/m2/s CO 26/211/252/2007 1 50
0 EPA16_CO__c3CO         $ROOT/NEI2016/v2021-06/2016fh_16j_emln_cmv_c3_12_0pt1degree_3D_month_$MM.ncf          CO  2002-2020/1-12/1/0 RF xyz  kg/m2/s CO 26/211/252/2007 1 50
0 EPA16_CO__pteguCO      $ROOT/NEI2016/v2021-06/2016fh_16j_emln_ptegu_0pt1degree_3D_month_$MM.ncf              CO  2002-2020/1-12/1/0 RF xyz  kg/m2/s CO 26/211/252/2007 1 50
0 EPA16_CO__ptogCO       $ROOT/NEI2016/v2021-06/2016fh_16j_emln_pt_oilgas_allinln_0pt1degree_3D_month_$MM.ncf  CO  2002-2020/1-12/1/0 RF xyz  kg/m2/s CO 26/211/252/2007 1 50
0 EPA16_CO__ptnonipmCO   $ROOT/NEI2016/v2021-06/2016fh_16j_emln_ptnonipm_allinln_0pt1degree_3D_month_$MM.ncf   CO  2002-2020/1-12/1/0 RF xyz  kg/m2/s CO 26/211/252/2007 1 50
0 EPA16_CO__onroad_canCO $ROOT/NEI2016/v2021-06/2016fh_16j_onroad_can_0pt1degree_month_$MM.ncf                 CO  2002-2020/1-12/1/0 RF xyz  kg/m2/s CO 26/211/252/2007 1 50
0 EPA16_CO__onroad_mexCO $ROOT/NEI2016/v2021-06/2016fh_16j_onroad_mex_0pt1degree_month_$MM.ncf                 CO  2002-2020/1-12/1/0 RF xyz  kg/m2/s CO 26/211/252/2007 1 50
0 EPA16_CO__otharCO      $ROOT/NEI2016/v2021-06/2016fh_16j_othar_0pt1degree_month_$MM.ncf                      CO  2002-2020/1-12/1/0 RF xyz  kg/m2/s CO 26/211/252/2007 1 50
0 EPA16_CO__othptCO      $ROOT/NEI2016/v2021-06/2016fh_16j_emln_othpt_0pt1degree_month_$MM.ncf                 CO  2002-2020/1-12/1/0 RF xyz  kg/m2/s CO 26/211/252/2007 1 50
)))NEI2016_MONMEAN

#==============================================================================
# --- CO: DICE-Africa emission inventory (Marais and Wiedinmyer, ES&T, 2016) ---
#
# DICE-Africa includes regional (Africa) emissions of biofuel and diffuse
# anthropogenic emissions from cars and motorcycles, biofuels, charcoal making
# and use, backup generators, agricultural waste burning for cooking, gas
# flares, and ad-hoc/informal oil refining.
#
# Other pollution sources (formal industry, power generation using fossil
# fuels) are from the EDGAR v4.3 inventory for CO, SO2, NH3, NOx BC, and OC.
#
# NMVOCs from sources not accounted for in DICE-Africa aren't included here,
# as these emissions are likely to be low compared to the DICE pollution
# sources and RETRO v1 as implemented in GEOS-Chem doesn't distinguish
# emissions by sector/activity.
#
# Emissions for 2013 are defined below, but DICE-Africa also includes
# emissions for 2006.  Developers recommend using population change to
# estimate emissions, if users want to use annual trends in pollutant
# emissions to estimate in other years.
#==============================================================================
(((DICE_Africa
0 DICE_CARS_CO          $ROOT/DICE_Africa/v2016-10/DICE-Africa-cars-2013-v01-4Oct2016.nc                       CO 2013/1/1/0 C xy g/m2/yr CO 26/2008     1 60
0 DICE_MOTORCYCLES_CO   $ROOT/DICE_Africa/v2016-10/DICE-Africa-motorcycles-2013-v01-4Oct2016.nc                CO 2013/1/1/0 C xy g/m2/yr CO 26/2008     1 60
0 DICE_BACKUPGEN_CO     $ROOT/DICE_Africa/v2016-10/DICE-Africa-generator-use-2013-v01-4Oct2016.nc              CO 2013/1/1/0 C xy g/m2/yr CO 26/2008     1 60
0 DICE_CHARCOALPROD_CO  $ROOT/DICE_Africa/v2016-10/DICE-Africa-charcoal-production-2013-v01-4Oct2016.nc        CO 2013/1/1/0 C xy g/m2/yr CO 26/2008/320 1 60
0 DICE_GASFLARE_CO      $ROOT/DICE_Africa/v2016-10/DICE-Africa-gas-flares-2013-v01-4Oct2016.nc                 CO 2013/1/1/0 C xy g/m2/yr CO 26/2008     1 60
0 DICE_AGBURNING_CO     $ROOT/DICE_Africa/v2016-10/DICE-Africa-household-crop-residue-use-2013-v01-4Oct2016.nc CO 2013/1/1/0 C xy g/m2/yr CO 26/2008     2 60
0 DICE_CHARCOALUSE_CO   $ROOT/DICE_Africa/v2016-10/DICE-Africa-charcoal-use-2013-v01-4Oct2016.nc               CO 2013/1/1/0 C xy g/m2/yr CO 26/2008     2 60
0 DICE_KEROSENE_CO      $ROOT/DICE_Africa/v2016-10/DICE-Africa-kerosene-use-2013-v01-4Oct2016.nc               CO 2013/1/1/0 C xy g/m2/yr CO 26/2008     1 60
0 DICE_OILREFINING_CO   $ROOT/DICE_Africa/v2016-10/DICE-Africa-adhoc-oil-refining-2006-v01-4Oct2016.nc         CO 2013/1/1/0 C xy g/m2/yr CO 26/2008     1 60
0 DICE_HOUSEFUELWOOD_CO $ROOT/DICE_Africa/v2016-10/DICE-Africa-household-fuelwood-use-2013-v01-4Oct2016.nc     CO 2013/1/1/0 C xy g/m2/yr CO 26/2008     2 60
0 DICE_OTHERFUELWOOD_CO $ROOT/DICE_Africa/v2016-10/DICE-Africa-other-fuelwood-use-2013-v01-4Oct2016.nc         CO 2013/1/1/0 C xy g/m2/yr CO 26/2008     2 60

# ---------------------------------------------------
#  Efficient Combustion Emissions from EDGAR
#  This makes up for sources that DICE-Africa lacks
# ---------------------------------------------------
0 AF_EDGAR_CO_POW $ROOT/EDGARv43/v2016-11/EDGAR_v43.CO.POW.0.1x0.1.nc emi_co  1970-2010/1/1/0 RF xy kg/m2/s CO 1201/26/52/2008 1 60
0 AF_EDGAR_CO_ENG $ROOT/EDGARv43/v2016-11/EDGAR_v43.CO.ENG.0.1x0.1.nc emi_co  1970-2010/1/1/0 RF xy kg/m2/s CO 1202/26/52/2008 1 60
0 AF_EDGAR_CO_IND $ROOT/EDGARv43/v2016-11/EDGAR_v43.CO.IND.0.1x0.1.nc emi_co  1970-2010/1/1/0 RF xy kg/m2/s CO 1203/26/52/2008 1 60
0 AF_EDGAR_CO_TNG $ROOT/EDGARv43/v2016-11/EDGAR_v43.CO.TNG.0.1x0.1.nc emi_co  1970-2010/1/1/0 RF xy kg/m2/s CO 1205/26/52/2008 1 60
0 AF_EDGAR_CO_PPA $ROOT/EDGARv43/v2016-11/EDGAR_v43.CO.PPA.0.1x0.1.nc emi_co  1970-2010/1/1/0 RF xy kg/m2/s CO 1207/26/52/2008 1 60
0 AF_EDGAR_CO_SWD $ROOT/EDGARv43/v2016-11/EDGAR_v43.CO.SWD.0.1x0.1.nc emi_co  1970-2010/1/1/0 RF xy kg/m2/s CO 1211/26/52/2008 1 60
)))DICE_Africa

#==============================================================================
# --- CO: CEDS v2 ---
#
# %%% This is the default global inventory. You may select either CEDS,
# EDGAR, HTAP or CMIP6_SFC_LAND_ANTHRO for the global base emissions %%%
#==============================================================================
(((CEDSv2
0 CEDS_CO_AGR  $ROOT/CEDS/v2021-06/$YYYY/CO-em-anthro_CMIP_CEDS_$YYYY.nc   CO_agr  1750-2019/1-12/1/0 C xy   kg/m2/s CO 26           1 5
0 CEDS_CO_ENE  $ROOT/CEDS/v2021-06/$YYYY/CO-em-anthro_CMIP_CEDS_$YYYY.nc   CO_ene  1750-2019/1-12/1/0 C xyL* kg/m2/s CO 26/37/35/315 1 5
0 CEDS_CO_IND  $ROOT/CEDS/v2021-06/$YYYY/CO-em-anthro_CMIP_CEDS_$YYYY.nc   CO_ind  1750-2019/1-12/1/0 C xyL* kg/m2/s CO 26/316       1 5
0 CEDS_CO_TRA  $ROOT/CEDS/v2021-06/$YYYY/CO-em-anthro_CMIP_CEDS_$YYYY.nc   CO_tra  1750-2019/1-12/1/0 C xy   kg/m2/s CO 26           1 5
0 CEDS_CO_RCO  $ROOT/CEDS/v2021-06/$YYYY/CO-em-anthro_CMIP_CEDS_$YYYY.nc   CO_rco  1750-2019/1-12/1/0 C xy   kg/m2/s CO 26           1 5
0 CEDS_CO_SLV  $ROOT/CEDS/v2021-06/$YYYY/CO-em-anthro_CMIP_CEDS_$YYYY.nc   CO_slv  1750-2019/1-12/1/0 C xy   kg/m2/s CO 26           1 5
0 CEDS_CO_WST  $ROOT/CEDS/v2021-06/$YYYY/CO-em-anthro_CMIP_CEDS_$YYYY.nc   CO_wst  1750-2019/1-12/1/0 C xy   kg/m2/s CO 26           1 5
)))CEDSv2

#==============================================================================
# --- CO: CEDS GBD-MAPS ---
#
# NOTES:
#  -- Reference: McDuffie et al. (2020, Earth System Science Data)
#  -- Anthropogenic source sectors: agriculture, energy, industry, road transport,
#       non-road/off-road transport, residential, commercial, other energy use,
#       solvents, waste, international shipping
#  -- Fuel categories: the combustion of total coal, solid biofuel, liquid oil
#       and gas, and all remaining sources
#       **To use, enable CEDS_byFuelType and CEDS_SHIP_byFuelType**
#==============================================================================
(((CEDS_GBDMAPS
>>>include $ROOT/CEDS/v2020-08/HEMCO_Config.CEDS_GBD-MAPS.rc
)))CEDS_GBDMAPS
(((CEDS_GBDMAPS_byFuelType
(((.not.CEDS_GBDMAPS.and..not.CEDSv2
>>>include $ROOT/CEDS/v2020-08/HEMCO_Config.CEDS_GBD-MAPS_byFuelType.rc
))).not.CEDS_GBDMAPS.and..not.CEDSv2
)))CEDS_GBDMAPS_byFuelType

#==============================================================================
# --- CO: CMIP6_SFC_LAND_ANTHRO ---
# CEDS (historical) or Shared Socioeconomic Pathways (future), consistent with
#   the CMIP6 simulation experimental design.
#
# Make sure that the desired $GCAP2SCENARIO is set above in SECTION SETTINGS
#==============================================================================
(((CMIP6_SFC_LAND_ANTHRO
0 CMIP6_CO_AGR  $ROOT/CMIP6/v2021-01/$GCAP2SCENARIO/$GCAP2SCENARIO_$YYYY.nc4 CO_agr 1750-2100/1-12/1/0 C xy kg/m2/s CO 26 1 5
0 CMIP6_CO_ENE  $ROOT/CMIP6/v2021-01/$GCAP2SCENARIO/$GCAP2SCENARIO_$YYYY.nc4 CO_ene 1750-2100/1-12/1/0 C xy kg/m2/s CO 26 1 5
0 CMIP6_CO_IND  $ROOT/CMIP6/v2021-01/$GCAP2SCENARIO/$GCAP2SCENARIO_$YYYY.nc4 CO_ind 1750-2100/1-12/1/0 C xy kg/m2/s CO 26 1 5
0 CMIP6_CO_TRA  $ROOT/CMIP6/v2021-01/$GCAP2SCENARIO/$GCAP2SCENARIO_$YYYY.nc4 CO_tra 1750-2100/1-12/1/0 C xy kg/m2/s CO 26 1 5
0 CMIP6_CO_RCO  $ROOT/CMIP6/v2021-01/$GCAP2SCENARIO/$GCAP2SCENARIO_$YYYY.nc4 CO_rco 1750-2100/1-12/1/0 C xy kg/m2/s CO 26 1 5
0 CMIP6_CO_SLV  $ROOT/CMIP6/v2021-01/$GCAP2SCENARIO/$GCAP2SCENARIO_$YYYY.nc4 CO_slv 1750-2100/1-12/1/0 C xy kg/m2/s CO 26 1 5
0 CMIP6_CO_WST  $ROOT/CMIP6/v2021-01/$GCAP2SCENARIO/$GCAP2SCENARIO_$YYYY.nc4 CO_wst 1750-2100/1-12/1/0 C xy kg/m2/s CO 26 1 5
)))CMIP6_SFC_LAND_ANTHRO

#==============================================================================
# CO: CEDS (historical) or Shared Socioeconomic Pathways (future) aircraft
# emissions, consistent with the CMIP6 simulation experimental design
#
# Make sure that the desired $GCAP2SCENARIO is set above in SECTION SETTINGS
#==============================================================================
(((CMIP6_AIRCRAFT
0 CMIP6_AIR_CO   $ROOT/CMIP6/v2021-01/$GCAP2SCENARIO/$GCAP2SCENARIO_$YYYY_AIR.$GCAP2VERTRESL.nc4 CO_air 1750-2100/1-12/1/0 C xyz kg/m2/s CO 26 20 1
)))CMIP6_AIRCRAFT

#==============================================================================
# --- CO: HTAP v2 ---
#
# %%% This is an optional inventory. You may select either CEDS, EDGAR,
#  or HTAP for the global base emissions %%%
#
# ==> HTAP ship emissions are listed in the ship emissions section below
# ==> Disable aircraft emissions and get them from AEIC instead.
#==============================================================================
(((HTAP
0 HTAP_CO_IND   $ROOT/HTAP/v2015-03/CO/EDGAR_HTAP_CO_INDUSTRY.generic.01x01.nc    emi_co 2008-2010/1-12/1/0 C xy kg/m2/s CO 506/528/26 1 4
0 HTAP_CO_POW   $ROOT/HTAP/v2015-03/CO/EDGAR_HTAP_CO_ENERGY.generic.01x01.nc      emi_co 2008-2010/1-12/1/0 C xy kg/m2/s CO 506/528/26 1 4
0 HTAP_CO_RES   $ROOT/HTAP/v2015-03/CO/EDGAR_HTAP_CO_RESIDENTIAL.generic.01x01.nc emi_co 2008-2010/1-12/1/0 C xy kg/m2/s CO 506/528/26 1 4
0 HTAP_CO_TRA   $ROOT/HTAP/v2015-03/CO/EDGAR_HTAP_CO_TRANSPORT.generic.01x01.nc   emi_co 2008-2010/1-12/1/0 C xy kg/m2/s CO 506/528/26 1 4
#0 HTAP_CO_AIR1  $ROOT/HTAP/v2015-03/CO/EDGAR_HTAP_CO_AIR_LTO.generic.01x01.nc     emi_co 2008-2010/1/1/0    C xy kg/m2/s CO 506/528/26 1 4
#0 HTAP_CO_AIR2  $ROOT/HTAP/v2015-03/CO/EDGAR_HTAP_CO_AIR_CDS.generic.01x01.nc     emi_co 2008-2010/1/1/0    C xy kg/m2/s CO 506/528/26 1 4
#0 HTAP_CO_AIR3  $ROOT/HTAP/v2015-03/CO/EDGAR_HTAP_CO_AIR_CRS.generic.01x01.nc     emi_co 2008-2010/1/1/0    C xy kg/m2/s CO 506/528/26 1 4
)))HTAP

#==============================================================================
# --- CO: Ship emissions ---
#
# ==> CEDS ship emissions are now the default.
# ==> If CEDS_SHIP is turned off above then ARCTAS should be used over ICOADS,
#     CORBETT, and HTAP for SO2 and ICOADS should be used for CO and NO.
# ==> Ship NO emissions are used by PARANOx and the extension number must be
#     adjusted accordingly. If PARANOx is turned off, set the ExtNr back to
#     zero.
#==============================================================================
(((SHIP
(((ICOADS_SHIP
0  ICOADS_SHIP_CO  $ROOT/ICOADS_SHIP/v2014-07/ICOADS.generic.1x1.nc CO 2002/1-12/1/0 C xy kg/m2/s  CO 6/10 10 2
)))ICOADS_SHIP

(((HTAP_SHIP
0 HTAP_SHIP_CO  $ROOT/HTAP/v2015-03/EDGAR_HTAP_CO_SHIPS.generic.01x01.nc CO 2008-2010/1/1/0 C xy kg/m2/s CO 506/528 10 4
)))HTAP_SHIP

(((CEDSv2_SHIP
0 CEDS_CO_SHP   $ROOT/CEDS/v2021-06/$YYYY/CO-em-anthro_CMIP_CEDS_$YYYY.nc CO_shp 1750-2019/1-12/1/0 C xy kg/m2/s CO 26 10 5
)))CEDSv2_SHIP

#==============================================================================
# --- CO: CMIP6_SHIP ---
# CEDS (historical) or Shared Socioeconomic Pathways (future), consistent with
# the CMIP6 simulation experimental design.
#
# Make sure that the desired $GCAP2SCENARIO is set above in SECTION SETTINGS
#==============================================================================
(((CMIP6_SHIP
0 CMIP6_CO_SHP  $ROOT/CMIP6/v2021-01/$GCAP2SCENARIO/$GCAP2SCENARIO_$YYYY.nc4 CO_shp 1750-2100/1-12/1/0 C xy kg/m2/s CO 26 10 5
)))CMIP6_SHIP

)))SHIP

#==============================================================================
# --- CO: AEIC 2019 aircraft emissions ---
#
# Data files are for 2019, but scale factors from 1990-2019 can be applied
# in order to get year-specific emissions.  See the notes in the AEIC2019
# scale factor section below for more information.
#==============================================================================
(((AEIC2019_DAILY
0 AEIC19_DAILY_CO    $ROOT/AEIC2019/v2022-03/2019/AEIC_2019$MM$DD.0.5x0.625.36L.nc CO 2019/1-12/1-31/0 C xyz kg/m2/s CO 241 20 1
)))AEIC2019_DAILY
(((AEIC2019_MONMEAN
0 AEIC19_MONMEAN_CO  $ROOT/AEIC2019/v2022-03/2019_monmean/AEIC_monmean_2019$MM.0.5x0.625.36L.nc CO 2019/1-12/1/0 C xyz kg/m2/s CO 241 20 1
)))AEIC2019_MONMEAN

#==============================================================================
# --- CO: RCP future emissions scenarios ---
#==============================================================================
(((RCP_3PD
0 RCP3PD_CO $ROOT/RCP/v2020-07/RCP_3PD/RCPs_anthro_CO_2005-2100_23474.nc ACCMIP 2005-2100/1/1/0 ID xy kg/m2/s CO - 1 1
)))RCP_3PD

(((RCP_45
0 RCP45_CO  $ROOT/RCP/v2020-07/RCP_45/RCPs_anthro_CO_2005-2100_27424.nc  ACCMIP 2005-2100/1/1/0 ID xy kg/m2/s CO - 1 1
)))RCP_45

(((RCP_60
0 RCP60_CO  $ROOT/RCP/v2020-07/RCP_60/RCPs_anthro_CO_2005-2100_43190.nc  ACCMIP 2005-2100/1/1/0 ID xy kg/m2/s CO - 1 1
)))RCP_60

(((RCP_85
0 RCP85_CO  $ROOT/RCP/v2020-07/RCP_85/RCPs_anthro_CO_2005-2100_43533.nc  ACCMIP 2005-2100/1/1/0 ID xy kg/m2/s CO - 1 1
)))RCP_85

#==============================================================================
# --- CO biomass burning (base emissions) ---
#==============================================================================

# --- QFED2 (v2.5r1) ---
(((QFED2
0 QFED_CO_PBL $ROOT/QFED/v2018-07/$YYYY/$MM/qfed2.emis_co.006.$YYYY$MM$DD.nc4 biomass 2000-2022/1-12/1-31/0/+12hour EF xyL=1:PBL     kg/m2/s CO 54/75/311 5 2
0 QFED_CO_FT  $ROOT/QFED/v2018-07/$YYYY/$MM/qfed2.emis_co.006.$YYYY$MM$DD.nc4 biomass 2000-2022/1-12/1-31/0/+12hour EF xyL=PBL:5500m kg/m2/s CO 54/75/312 5 2
)))QFED2

# --- GFAS ---
(((GFAS
0 GFAS_CO $ROOT/GFAS/v2018-09/$YYYY/GFAS_$YYYY$MM.nc cofire 2003-2021/1-12/1-31/0 C xyL=1:scal300 kg/m2/s CO 75 5 3
)))GFAS

# --- BB4MIPs ---
(((BB4MIPS
0 CMIP6_BB_CO $ROOT/CMIP6/v2021-01/$GCAP2SCENARIO/$GCAP2SCENARIO_$YYYY.nc4 CO_bbn 1750-2100/1-12/1/0 C xyL=1:PBL kg/m2/s CO 75 5 3
)))BB4MIPS

)))USE_CO_DATA

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%                        CO2 EMISSIONS                               %%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
(((USE_CO2_DATA

#==============================================================================
# --- CO2: FOSSIL FUEL EMISSIONS ---
#
# National CO2 fossil fuel emissions exclude international shipping and aviation
# Carbon Dioxide Information Analysis Center (CDIAC) 1x1 national emissions
# (Andres et al., 2011)
# Open-source Data Inventory of Anthropogenic CO2 (ODIAC) 1x1 national emissions
# (Oda & Maksyutov, 2011)
#
# ---> Recommended option: ODIAC (set FOSSIL_ODIAC = true)
#      ODIAC updated to v2019 by J. Fisher and Y. Cao, 12/2019
#==============================================================================
(((FOSSIL_CDIAC
0 FOSSILCO2_CDIAC   $ROOT/CO2/v2014-09/FOSSIL/CDIAC_v2014.monthly.generic.1x1.nc CO2 1980-2014/1-12/1/0 C xy kg/m2/s CO2 440/441/480 1 1
)))FOSSIL_CDIAC

(((FOSSIL_ODIAC
0 FOSSILCO2_ODIAC   $ROOT/CO2/v2022-11/FOSSIL/ODIAC_CO2.monthly.generic.1x1.nc   CO2 2000-2018/1-12/1/0 C xy kg/m2/s CO2 440/441/480 1 2
)))FOSSIL_ODIAC

#==============================================================================
# --- CO2: OCEAN EXCHANGE EMISSIONS ---
#
# Ocean uptake/emission from Takahashi et al. (2009) and earlier works.
# Climatological exchange for 2000 (-1.4 PgC/yr) or scaled ocean exchange, which
# accounts for changing atmospheric CO2 (as in the paper) are now available.
#
# ---> Recommended option: scaled ocean exchange (set OCEAN_EXCH_SCALED = true)
#==============================================================================
(((OCEAN_EXCH_TAKA09
0 OCEANCO2_TAKA_ANNUAL    $ROOT/CO2/v2022-11/OCEAN/Taka2009_CO2_Annual.nc      CO2 2000/1/1/0    C xy kg/m2/s CO2   - 2 1
0 OCEANCO2_TAKA_MONTHLY   $ROOT/CO2/v2022-11/OCEAN/Taka2009_CO2_Monthly.nc     CO2 2000/1-12/1/0 C xy kg/m2/s CO2   - 2 2
)))OCEAN_EXCH_TAKA09

(((OCEAN_EXCH_SCALED
0 OCEANCO2_SCALED_MONTHLY  $ROOT/CO2/v2022-11/OCEAN/Scaled_Ocean_CO2_monthly.nc CO2 2000-2013/1-12/1/0 C xy kg/m2/s CO2   - 2 3
)))OCEAN_EXCH_SCALED

#==============================================================================
# --- CO2: BALANCED BIOSPHERE EXCHANGE ---
#
# These emissions have a seasonal cycle of uptake/emission but have a net
# annual uptake of close to zero (balanced).
# CASA model daily emissions for a single year at coarse resolution have a
# diurnal cycle imposed (Olsen & Randerson, 2004)
# SiB3 emissions are from Nick Parazoo (Messerschmidt et al., 2011) for
# 2006-2010.
#
# ---> Recommended option: SiB3 (set BBIO_SIB3 = true)
#==============================================================================
(((BBIO_DIURNAL
0 BBIOCO2_DIURNAL $ROOT/CO2/v2014-09/BBIO/BBIO_diurnal_CO2.nc CO2 1985/1-12/1-31/0-23 C xy kg/m2/s CO2 - 3 1
)))BBIO_DIURNAL

(((BBIO_SIB3
0 SIB_BBIO_CO2    $ROOT/CO2/v2022-11/BIO/SiB3_3hr_NEP.nc CO2 2006-2010/1-12/1-31/0-23 C xy kg/m2/s CO2 - 3 1
)))BBIO_SIB3

#==============================================================================
# --- CO2: NET TERRESTRIAL EXCHANGE ---
#
# TransCom annual net/residual terrestrial biosperhic CO2 (Baker et al., 2006)
#
# ---> Recommended for use in forward modelling, optional for
#      inversion/assimilation
#==============================================================================
(((NET_TERR_EXCH
0 CO2_NET_TERRESTRIAL  $ROOT/CO2/v2022-11/BIO/Net_terrestrial_exch_5.29Pg.generic.1x1.nc CO2 2000/1/1/0 C xy kg/m2/s CO2   - 5 1
)))NET_TERR_EXCH

#==============================================================================
# --- CO2: SHIP EMISSIONS ---
#==============================================================================
(((SHIP

(((ICOADS_SHIP
# The spatial distribution is scaled with global annual scale factors #50
0 ICOADS_CO2_SHIP    $ROOT/ICOADS_SHIP/v2014-07/ICOADS_ship_CO2_2004.generic.1x1.nc CO2 2004/1-12/1/0 C xy kg/m2/s CO2 50 6 1
)))ICOADS_SHIP

(((CEDSv2_SHIP
0 CEDS_CO2_SHP  $ROOT/CEDS/v2018-08/$YYYY/CO2-em-anthro_CMIP_CEDS_$YYYY.nc  CO2_shp  1750-2014/1-12/1/0 C xy kg/m2/s CO2 - 6 1
)))CEDSv2_SHIP

)))SHIP

#==============================================================================
# --- CO2: AEIC 2019 aircraft emissions ---
#
# Data files are for 2019, but scale factors from 1990-2019 can be applied
# in order to get year-specific emissions.  See the notes in the AEIC2019
# scale factor section below for more information.
#==============================================================================
(((AEIC2019_DAILY
0 AEIC19_DAILY_CO2  $ROOT/AEIC2019/v2022-03/2019/AEIC_2019$MM$DD.0.5x0.625.36L.nc FUELBURN 2019/1-12/1-31/0 C xyz kg/m2/s CO2 241/260 20 1
)))AEIC2019_DAILY

(((AEIC2019_MONMEAN
0 AEIC19_MONMEAN_CO2  $ROOT/AEIC2019/v2022-03/2019_monmean/AEIC_monmean_2019$MM.0.5x0.625.36L.nc FUELBURN 2019/1-12/1/0 C xyz kg/m2/s CO2 241/260 20 1
)))AEIC2019_MONMEAN

#==============================================================================
# --- CO2: SURFACE CORRECTION FOR CO OXIDATION ---
#
# These emissions will be subtracted!
# Fossil fuel CO and CH4 are based on the national distribution scaled with #10.
# Biogenic CH4 data is converted to CO2 using scale factor #20.
# Isoprene and monoterpene are converted from kgC to kgCO2 by scale factor #21.
# Approach is described in Nassar et al. (2010)
#
# ---> Recommended for use if chemical production (above) is turned on
#
# Now multiply by scale factor of -1.0 (#1 in scale factor section below)
# to make these emissions negative, so that they will be subtracted.
===============================================================================
(((CO2CORR
0 FOSSILCO2_MONTHLY $ROOT/CO2/v2022-11/FOSSIL/ODIAC_CO2.monthly.generic.1x1.nc CO2     2000-2018/1-12/1/0 C xy kg/m2/s CO2 410/440/441/480/1 8 1
0 CO2_LIVESTOCK     $ROOT/CO2/v2022-11/CHEM/CH4_source.geos.2x25.nc            CH4_004 2004/1-12/1/0      C xy kg/m2/s CO2 420/1             8 1
0 CO2_WASTE         $ROOT/CO2/v2022-11/CHEM/CH4_source.geos.2x25.nc            CH4_005 2004/1-12/1/0      C xy kg/m2/s CO2 420/1             8 1
0 CO2_RICE          $ROOT/CO2/v2022-11/CHEM/CH4_source.geos.2x25.nc            CH4_007 2004/1-12/1/0      C xy kg/m2/s CO2 420/1             8 1
0 CO2_WETLANDS      $ROOT/CO2/v2022-11/CHEM/CH4_source.geos.2x25.nc            CH4_010 2004/1-12/1/0      C xy kg/m2/s CO2 420/1             8 1
0 CO2_NATURAL       $ROOT/CO2/v2022-11/CHEM/CH4_source.geos.2x25.nc            CH4_012 2004/1-12/1/0      C xy kg/m2/s CO2 420/1             8 1
0 CO2_ISOPRENE      $ROOT/CO2/v2022-11/CHEM/Isoprene-2004.geos.2x25.nc         ISOP    2004/1-12/1/0      C xy kg/m2/s CO2 421/430/1         8 1
0 CO2_MONOTERP      $ROOT/CO2/v2022-11/CHEM/Monoterpene-2004.geos.2x25.nc      MONOT   2004/1-12/1/0      C xy kg/m2/s CO2 421/430/1         8 1
)))CO2CORR

)))USE_CO2_DATA

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%                OCS FLUXES (cf Kevin Bowman/JPL)                    %%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

(((USE_OCS_DATA

(((OCS_ANTHRO_FLUX
0 OCS_ANTHRO $ROOT/OCS_FLUX/v2022-11/anthro_v3/v3_anthro/$YYYY/$MM.nc COS_Flux 2012-2016/1-12/1/0 C xy kgS/km2/s OCS 70/71 1 1
)))OCS_ANTHRO_FLUX

(((OCS_BIOMASS_FLUX
0 OCS_BIOMASS $ROOT/OCS_FLUX/v2022-11/Biomassburn/CMSV2-CO2-GFED3-ZC/$YYYY/$MM.nc COS_Flux 2012-2016/1-12/1/0 C xy kgS/km2/s OCS 70/71 2 1
)))OCS_BIOMASS_FLUX

(((OCS_MISSING_OCEAN_FLUX
0 OCS_MISSING_OCEAN $ROOT/OCS_FLUX/v2022-11/MissingOcean/LUKAI/$YYYY/$MM.nc COS_Flux 2012-2016/1-12/1/0 C xy kgS/km2/s OCS 70/71 3 1
)))OCS_MISSING_OCEAN_FLUX

(((OCS_OCEAN_FLUX
0 OCS_OCEAN $ROOT/OCS_FLUX/v2022-11/OceanCOS/Kettle/$YYYY/$MM.nc COS_Flux 2012-2016/1-12/1/0 C xy kgS/km2/s OCS 70/71 4 1
)))OCS_OCEAN_FLUX

)))USE_OCS_DATA

###############################################################################
### EXTENSION DATA (subsection of BASE EMISSIONS SECTION)
###
### These fields are needed by the extensions listed above. The assigned ExtNr
### must match the ExtNr entry in section 'Extension switches'. These fields
### are only read if the extension is enabled.  The fields are imported by the
### extensions by field name.  The name given here must match the name used
### in the extension's source code.
###############################################################################

#==============================================================================
# --- GFED biomass burning emissions (Extension 111)
# NOTE: These are the base emissions in kgDM/m2/s.
#==============================================================================
(((USE_CH4_DATA.or.USE_CO2_DATA

(((GFED4
111 GFED_TEMP       $ROOT/GFED4/v2020-02/$YYYY/GFED4_gen.025x025.$YYYY$MM.nc  DM_TEMP       1997-2019/1-12/01/0    RF xy kgDM/m2/s * - 1 1
111 GFED_AGRI       $ROOT/GFED4/v2020-02/$YYYY/GFED4_gen.025x025.$YYYY$MM.nc  DM_AGRI       1997-2019/1-12/01/0    RF xy kgDM/m2/s * - 1 1
111 GFED_DEFO       $ROOT/GFED4/v2020-02/$YYYY/GFED4_gen.025x025.$YYYY$MM.nc  DM_DEFO       1997-2019/1-12/01/0    RF xy kgDM/m2/s * - 1 1
111 GFED_BORF       $ROOT/GFED4/v2020-02/$YYYY/GFED4_gen.025x025.$YYYY$MM.nc  DM_BORF       1997-2019/1-12/01/0    RF xy kgDM/m2/s * - 1 1
111 GFED_PEAT       $ROOT/GFED4/v2020-02/$YYYY/GFED4_gen.025x025.$YYYY$MM.nc  DM_PEAT       1997-2019/1-12/01/0    RF xy kgDM/m2/s * - 1 1
111 GFED_SAVA       $ROOT/GFED4/v2020-02/$YYYY/GFED4_gen.025x025.$YYYY$MM.nc  DM_SAVA       1997-2019/1-12/01/0    RF xy kgDM/m2/s * - 1 1

(((GFED_daily
111 GFED_FRAC_DAY   $ROOT/GFED4/v2020-02/$YYYY/GFED4_dailyfrac_gen.025x025.$YYYY$MM.nc GFED_FRACDAY 2003-2019/1-12/1-31/0  RF xy 1 * - 1 1
)))GFED_daily

(((GFED_3hourly
111 GFED_FRAC_3HOUR $ROOT/GFED4/v2020-02/$YYYY/GFED4_3hrfrac_gen.025x025.$YYYY$MM.nc   GFED_FRAC3HR 2003-2019/1-12/1/0-23  RF xy 1 * - 1 1
)))GFED_3hourly
)))GFED4

)))USE_CH4_DATA.or.USE_CO2_DATA

)))EMISSIONS

###############################################################################
### NON-EMISSIONS DATA (subsection of BASE EMISSIONS SECTION)
###
### Non-emissions data. The following fields are read through HEMCO but do
### not contain emissions data. The extension number is set to wildcard
### character denoting that these fields will not be considered for emission
### calculation. A given entry is only read if the assigned species name is
### an HEMCO species.
###############################################################################

#==============================================================================
# --- Time zones (offset to UTC) ---
#==============================================================================
* TIMEZONES $ROOT/TIMEZONES/v2015-02/timezones_voronoi_1x1.nc UTC_OFFSET 2000/1/1/0 C xy count * - 1 1

#==============================================================================
# --- Meteorology fields ---
#==============================================================================
(((METEOROLOGY

>>>include ${RUNDIR_MET_FIELD_CONFIG}

)))METEOROLOGY

#==============================================================================
# --- GEOS-Chem restart file ---
#==============================================================================
(((GC_RESTART
* SPC_           ./Restarts/GEOSChem.Restart.$YYYY$MM$DD_$HH$MNz.nc4 SpeciesRst_?ALL?    $YYYY/$MM/$DD/$HH EFYO xyz 1 * - 1 1
* DELPDRY        ./Restarts/GEOSChem.Restart.$YYYY$MM$DD_$HH$MNz.nc4 Met_DELPDRY         $YYYY/$MM/$DD/$HH EY   xyz 1 * - 1 1
)))GC_RESTART

#==============================================================================
# --- GEOS-Chem boundary condition file ---
#==============================================================================
(((GC_BCs
* BC_ $ROOT/SAMPLE_BCs/v2021-07/CH4/GEOSChem.BoundaryConditions.$YYYY$MM$DD_0000z.nc4 SpeciesBC_?ADV?  1980-2021/1-12/1-31/* EFY xyz 1 * - 1 1
)))GC_BCs

#==============================================================================
# --- Olson land map masks ---
#==============================================================================
(((OLSON_LANDMAP
* LANDTYPE00 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE00 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE01 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE01 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE02 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE02 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE03 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE03 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE04 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE04 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE05 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE05 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE06 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE06 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE07 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE07 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE08 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE08 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE09 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE09 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE10 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE10 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE11 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE11 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE12 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE12 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE13 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE13 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE14 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE14 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE15 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE15 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE16 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE16 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE17 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE17 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE18 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE18 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE19 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE19 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE20 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE20 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE21 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE21 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE22 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE22 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE23 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE23 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE24 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE24 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE25 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE25 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE26 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE26 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE27 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE27 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE28 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE28 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE29 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE29 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE30 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE30 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE31 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE31 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE32 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE32 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE33 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE33 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE34 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE34 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE35 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE35 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE36 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE36 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE37 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE37 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE38 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE38 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE39 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE39 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE40 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE40 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE41 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE41 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE42 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE42 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE43 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE43 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE44 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE44 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE45 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE45 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE46 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE46 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE47 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE47 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE48 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE48 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE49 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE49 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE50 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE50 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE51 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE51 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE52 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE52 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE53 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE53 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE54 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE54 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE55 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE55 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE56 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE56 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE57 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE57 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE58 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE58 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE59 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE59 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE60 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE60 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE61 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE61 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE62 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE62 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE63 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE63 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE64 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE64 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE65 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE65 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE66 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE66 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE67 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE67 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE68 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE68 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE69 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE69 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE70 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE70 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE71 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE71 1985/1/1/0 C xy 1 * - 1 1
* LANDTYPE72 $ROOT/OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc LANDTYPE72 1985/1/1/0 C xy 1 * - 1 1
)))OLSON_LANDMAP

#==============================================================================
# --- Yuan processed MODIS leaf area index data ---
#
# Source: Yuan et al 2011, doi:10.1016/j.rse.2011.01.001
#         http://globalchange.bnu.edu.cn/research/lai
#
# NOTES:
# (1) LAI data corresponding to each Olson land type is stored in
#      separate netCDF variables (XLAI00, XLAI01, ... XLAI72).
#      The "XLAI" denotes that the files are prepared in this way.
# (2) Units are "cm2 leaf/cm2 grid box".
# (3) Data is timestamped every 8 days, starting from the 2nd of the month.
#==============================================================================
(((YUAN_MODIS_LAI
* XLAI00 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI00 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI01 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI01 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI02 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI02 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI03 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI03 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI04 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI04 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI05 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI05 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI06 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI06 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI07 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI07 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI08 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI08 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI09 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI09 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI10 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI10 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI11 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI11 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI12 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI12 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI13 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI13 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI14 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI14 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI15 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI15 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI16 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI16 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI17 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI17 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI18 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI18 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI19 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI19 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI20 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI20 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI21 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI21 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI22 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI22 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI23 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI23 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI24 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI24 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI25 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI25 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI26 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI26 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI27 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI27 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI28 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI28 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI29 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI29 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI30 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI30 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI31 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI31 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI32 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI32 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI33 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI33 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI34 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI34 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI35 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI35 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI36 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI36 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI37 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI37 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI38 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI38 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI39 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI39 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI40 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI40 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI41 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI41 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI42 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI42 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI43 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI43 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI44 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI44 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI45 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI45 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI46 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI46 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI47 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI47 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI48 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI48 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI49 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI49 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI50 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI50 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI51 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI51 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI52 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI52 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI53 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI53 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI54 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI54 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI55 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI55 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI56 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI56 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI57 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI57 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI58 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI58 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI59 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI59 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI60 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI60 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI61 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI61 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI62 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI62 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI63 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI63 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI64 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI64 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI65 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI65 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI66 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI66 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI67 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI67 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI68 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI68 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI69 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI69 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI70 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI70 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI71 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI71 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
* XLAI72 $ROOT/Yuan_XLAI/v2021-06/Yuan_proc_MODIS_XLAI.025x025.$YYYY.nc XLAI72 2000-2020/1-12/1-31/0 I xy cm2/cm2 * - 1 1
)))YUAN_MODIS_LAI

#==============================================================================
# --- Quantities needed for the analytical inversion (via IMI) ---
#==============================================================================
(((USE_CH4_DATA

#------------------------------------------------------------------------------
# --- Files needed for analytical inversion ---
#
# These fields are are only used if analytical_inversion?' is activated in
# geoschem_config.yml. These fields are obtained from HEMCO and applied in
# GEOS-Chem/GeosCore/global_ch4_mod.F90.
#
# Entries below provided for examples only. Add your own here!
#------------------------------------------------------------------------------
(((AnalyticalInv

# State vector file
* CH4_STATE_VECTOR StateVector.nc StateVector 2009/1/1/0 C xy 1 * - 1 1

)))AnalyticalInv

#------------------------------------------------------------------------------
# --- Scale factors for posterior run ---
#
# Enable emission scale factors by setting the use_emission_scale_factor or
# use_OH_scale_factor options to true in geoschem_config.yml. These fields are
# obtained from HEMCO and applied in GEOS-Chem/GeosCore/global_ch4_mod.F90.
#
# Entries below are provided for examples only. Add your own here!
#------------------------------------------------------------------------------
(((Emis_ScaleFactor
* EMIS_SF    gridded_posterior.nc  ScaleFactor  2000/1/1/0       C xy 1 * - 1 1
)))Emis_ScaleFactor

(((OH_ScaleFactor
* OH_SF      Post_SF_OH.nc         SF_OH        2010-2017/1/1/0  E xy 1 * - 1 1
)))OH_ScaleFactor

)))USE_CH4_DATA

#==============================================================================
# --- Chemistry inputs ---
#==============================================================================
(((CHEMISTRY_INPUT

#------------------------------------------------------------------------------
# --- Global OH fields ---
#
# If GLOBAL_OH_GCv5 is turned on:
#   - Use GEOS-Chem v5-07 OH if GLOBAL_OH_GCv5 is turned on.
#   - This is the recommended setting for CH4 and the IMI.
#
# If GLOBAL_OH_GCv5 is turned off:
#   - Use the OH archived from the last 10-yr benchmark simulation.
#------------------------------------------------------------------------------
(((GLOBAL_OH

# --- OH from GEOS-Chem v5-07 [kg/m3] ---
(((GLOBAL_OH_GCv5
* GLOBAL_OH  $ROOT/OH/v2014-09/v5-07-08/OH_3Dglobal.geos5.47L.4x5.nc OH  1985/1-12/1/0 C xyz kg/m3 * - 1 1
)))GLOBAL_OH_GCv5

# --- OH from the last 10-yr benchmark [mol/mol dry air] ---
(((.not.GLOBAL_OH_GCv5
${RUNDIR_GLOBAL_OH}
))).not.GLOBAL_OH_GCv5

)))GLOBAL_OH

#------------------------------------------------------------------------------
# --- Quantities needed for CH4 chemistry ---
#------------------------------------------------------------------------------
(((USE_CH4_DATA

# --- Global CH4 loss, Cl fields for CH4 simulation [1/s] -----------
(((GMI_CH4_LOSS
${RUNDIR_CH4_LOSS}
)))GMI_CH4_LOSS

# --- Cl from the last 10-yr benchmark [mol/mol dry air] ---
(((GLOBAL_CL
${RUNDIR_GLOBAL_Cl}
)))GLOBAL_CL

)))USE_CH4_DATA

#------------------------------------------------------------------------------
# --- Quantities needed for CO chemistry ---
#------------------------------------------------------------------------------
(((USE_CO_DATA

# -- P(CO) from CH4 and NMVOC from the last 10-yr benchmark [molec/cm3/s] ---
(((PROD_CO_CH4
${RUNDIR_PCO_CH4}
)))PROD_CO_CH4
(((PROD_CO_NMVOC
${RUNDIR_PCO_NMVOC}
)))PROD_CO_NMVOC

# --- GMI chemistry: prod/loss rates (for strato-/mesosphere) ---
# --- Units: prod [v/v/s]; loss [1/s]                         ---
(((GMI_PROD_LOSS
${RUNDIR_GMI_LOSS_CO}
${RUNDIR_GMI_PROD_CO}
)))GMI_PROD_LOSS

)))USE_CO_DATA

#------------------------------------------------------------------------------
# --- Quantities needed for CO2 chemistry ---
#------------------------------------------------------------------------------

# --- CHEMICAL PRODUCTION FROM CO OXIDATION ------
# --- Recommended for use in forward modelling ---
# --- Optional for inversion/assimilation --------
(((USE_CO2_DATA

(((CO2_COPROD
${RUNDIR_CO2_COPROD}
)))CO2_COPROD

)))USE_CO2_DATA

)))CHEMISTRY_INPUT

### END SECTION BASE EMISSIONS ###

###############################################################################
### BEGIN SECTION SCALE FACTORS
###############################################################################

# ScalID Name sourceFile sourceVar sourceTime C/R/E SrcDim SrcUnit Oper

(((EMISSIONS

#==============================================================================
# --- Multiply by -1 to get a "negative" flux.
#==============================================================================
1 NEGATIVE -1.0 - - - xy 1 1

#==============================================================================
# --- Seasonal scaling factors ----
#==============================================================================
(((EDGARv6.or.GEPA
10 MANURE_SF $ROOT/CH4/v2017-10/Seasonal_SF/EMICH4_Manure_ScalingFactors.WithClimatology.nc  sf_ch4 2008-2016/1-12/1/0 C xy 1 1
11 RICE_SF   $ROOT/CH4/v2017-10/Seasonal_SF/EMICH4_Rice_ScalingFactors.SetMissing.nc         sf_ch4 2012/1-12/1/0  C xy 1 1
)))EDGARv6.or.GEPA

#==============================================================================
# --- Diurnal scale factors ---
#==============================================================================
26 GEIA_TOD_FOSSIL 0.45/0.45/0.6/0.6/0.6/0.6/1.45/1.45/1.45/1.45/1.4/1.4/1.4/1.4/1.45/1.45/1.45/1.45/0.65/0.65/0.65/0.65/0.45/0.45 - - - xy 1 1

#==============================================================================
# --- Seasonal scale factors ---
#==============================================================================
(((DICE_Africa
# from GEIA:
30 GEIA_SEASON_NOX  $ROOT/GEIA/v2014-07/GEIA_monthscal.generic.1x1.nc NOXrat 1985/1-12/1/0 C xy unitless 1
31 GEIA_SEASON_SO2  $ROOT/GEIA/v2014-07/GEIA_monthscal.generic.1x1.nc SO2rat 1985/1-12/1/0 C xy unitless 1
)))DICE_Africa

#==============================================================================
# Diurnal variablity of Chinese power plants from Liu et al. (EST, 2019)
#==============================================================================
33 PKU_pow_NOx              0.94/0.93/0.93/0.92/0.93/0.96/0.98/0.99/0.98/1.00/1.05/1.03/1.01/1.04/1.05/1.05/1.07/1.06/1.05/1.05/1.03/1.01/0.99/0.96 - - - xy unitless 1 1009
34 PKU_pow_SO2              0.92/0.91/0.90/0.89/0.90/0.93/0.97/0.98/0.99/1.02/1.07/1.04/1.03/1.06/1.06/1.07/1.10/1.10/1.06/1.05/1.03/1.01/0.99/0.93 - - - xy unitless 1 1009
35 PKU_pow_PM_BC_POC_VOC_CO 0.95/0.93/0.92/0.91/0.90/0.93/0.97/0.97/0.99/1.03/1.04/1.03/1.02/1.03/1.05/1.07/1.07/1.07/1.06/1.04/1.03/1.02/1.00/0.97 - - - xy unitless 1 2009

### These scale factors undo (Oper=-1) the global diurnal scale factors over China (Mask=1009)
#36 EDGAR_TODNOX_UNDO $ROOT/EDGARv42/v2015-02/NO/EDGAR_hourly_NOxScal.nc NOXscale 2000/1/1/* C xy unitless -1 2009
37 GEIA_TOD_FOSSIL_UNDO 0.45/0.45/0.6/0.6/0.6/0.6/1.45/1.45/1.45/1.45/1.4/1.4/1.4/1.4/1.45/1.45/1.45/1.45/0.65/0.65/0.65/0.65/0.45/0.45 - - - xy unitless -1 2009

#==============================================================================
# --- Scale factors used for species conversions ---
#==============================================================================

# VOC speciations
(((APEI.or.EDGARv43.or.DICE_Africa.or.QFED2
52 COPROD_FOSSIL    1.02  - - - xy unitless 1
54 COPROD_BIOMASS   1.05  - - - xy unitless 1
)))APEI.or.EDGARv43.or.DICE_Africa.or.QFED2

# S to OCS ratio
70 StoOCS  MATH:60.07/32.06 - - - xy unitless 1

# km-2 to m-2
71 Km2ToM2 1.0e-6 - - - xy unitless 1

#==============================================================================
# --- QFED2 diurnal scale factors ---
#
# Fire diurnal scaling factors (% per hour) from WRAP Report to Project
# No. 178-6, July 2005
#==============================================================================
(((QFED2.or.GFAS.or.BB4MIPS
75 QFED2_TOD 0.1392/0.1392/0.1368/0.1368/0.1368/0.1368/0.1368/0.1368/0.1368/0.48/0.96/1.68/2.4/3.12/3.84/4.08/2.88/1.68/0.96/0.1368/0.1368/0.1368/0.1368/0.1368 - - - xy unitless 1
)))QFED2.or.GFAS.or.BB4MIPS

#==============================================================================
# --- EPA NEI day-of-week scale factors ---
#==============================================================================
(((NEI2016_MONMEAN
211 NEI99_DOW_CO  $ROOT/NEI2005/v2023-02/NEI99.dow.geos.corrected.012023.1x1.nc CO 1999/1-12/WD/0 C xy 1 1
)))NEI2016_MONMEAN
#==============================================================================
# --- EPA NEI2016 annual scale factors ---
#
# Annual scale factors were computed from the EPA Trends Report for Tier 1 CAPS
# (obtained 21 Sep 2021) using the "Total without wildfires" field. NH3 and
# PM2.5 only had values for 2002-2020, while the remaining species had yearly
# values for 1990-2020 (we include only 2002-2020 here).
#
# See NEI2016/v2021-06/national_tier1_caps+HEMCOscaling.xlsx for details.
#==============================================================================
(((NEI2016_MONMEAN
252 NEI2016_CO_YRSCALE 1.817/1.767/1.716/1.666/1.610/1.554/1.393/1.251/1.272/1.261/1.220/1.178/1.137/1.095/1.000/0.973/0.950/0.927/0.904 - 2002-2020/1/1/0 C xy 1 1
)))NEI2016_MONMEAN

#==============================================================================
# --- AEIC2019 aircraft emissions scale factors ---
#
# See http://geoschemdata.wustl.edu/ExtData/HEMCO/AEIC2019/v2022-03/AEIC_2019_technical_note.pdf
#==============================================================================
(((AEIC2019_DAILY.or.AEIC2019_MONMEAN

#------------------------------------------------------------------------------
# Scaling factors for 1990-2019 derived from Lee et al. (2021).  Lee et al.
# (2021) only covers 1990 to 2018, so to get to 2019 it is assumed that the
# growth from 2017 to 2018 is the same as that from 2018 to 2019.
# So the formula is something like:
#
# Emissions of CO in 2009 = AEIC 2019 emissions of CO
#                          * (Lee 2017 CO        / Lee 2018 fuel burn)
#                          * (Lee 2009 fuel burn / Lee 2018 fuel burn)
#
# So in this case, we use the Lee 2017/Lee 2018 value to scale AEIC’s
# emissions to the “2018” values, and then scale directly using the Lee et al
# fuel burn. This ensures that, when running with year 2019, you get an
# unadjusted version of the AEIC2019 inventory, and all previous years are
# scaled down.
#
# All scaling factors are included in here in HEMCO_Config.rc.
#------------------------------------------------------------------------------
(((AEIC_SCALE_1990_2019
241 AC_FBMULT  0.506/0.489/0.490/0.493/0.517/0.529/0.553/0.570/0.581/0.600/0.631/0.607/0.608/0.608/0.646/0.678/0.686/0.706/0.703/0.666/0.700/0.721/0.728/0.749/0.773/0.815/0.854/0.905/0.952/1.000 - 1990-2019/1/1/0 C xy 1 1
)))AEIC_SCALE_1990_2019

# If not applying 1990-2019 scale factors, use 1.0
(((.not.AEIC_SCALE_1990_2019
241 AC_FBMULT 1.000000e+0 - -  - xy 1 1
))).not.AEIC_SCALE_1990_2019

#------------------------------------------------------------------------------
# Assume 3.159 kg CO2 from every kg of fuel burned
# cf Hileman, Stratton, & Donohoo, _J. Propul. Power_, 26(6), 1184–1196, 2010.
#------------------------------------------------------------------------------
260 AEIC19_FBtoCO2 3.159 - - - xy unitless 1

)))AEIC2019_DAILY.or.AEIC2019_MONMEAN

#==============================================================================
# --- GFAS scale factors ---
#==============================================================================
(((GFAS
300 GFAS_EMITL $ROOT/GFAS/v2018-09/$YYYY/GFAS_$YYYY$MM.nc mami 2003-2021/1-12/1-31/0 C xy m 1
)))GFAS

#==============================================================================
# --- QFED vertical partitioning ---
# Following Fischer et al. (2014) and Travis et al. (2016), emit 35% of QFED
# emissions above the PBL.
#==============================================================================
(((QFED2
311 QFED_PBL_FRAC 0.65 - - - xy 1 1
312 QFED_FT_FRAC  0.35 - - - xy 1 1
)))QFED2

#==============================================================================
# --- CEDS vertical partitioning ---
#==============================================================================
(((CEDSv2.or.CEDS_GBDMAPS
315 ENERGY_LEVS   $ROOT/VerticalScaleFactors/v2021-05/gc_layers.nc g_energy   2017/1/1/0 C xyz 1 1
316 INDUSTRY_LEVS $ROOT/VerticalScaleFactors/v2021-05/gc_layers.nc g_industry 2017/1/1/0 C xyz 1 1
317 SHIP_LEVS     $ROOT/VerticalScaleFactors/v2021-05/gc_layers.nc cmv_c3     2017/1/1/0 C xyz 1 1
)))CEDSv2.or.CEDS_GBDMAPS

#==============================================================================
# --- DICE-Africa ---
#==============================================================================
(((DICE_Africa
# Charcoal production scale factor to reduce charcoal production
# by a factor of 5 after finding error in implementation of emission factors.
320 DICE_CP_SF    0.20 - - - xy 1 1
)))DICE_Africa

#==============================================================================
# --- FOSSIL FUEL CO2 SURFACE CORRECTION FACTOR ---
#
# Fossil fuel CO2 emissions must be scaled down to avoid counting emissions
# already accounted for by CO and CH4 oxidation to CO2.
#
# The necssary annual reduction in emissions related to fossil fuels is
# determined as:  global C mass from FF (CO+CH4) / global C mass from FF CO2
#
# Note: GEOS-Chem v8-03-02 to v9-02, assumed: 0.0489 for all years but that
# value seems to correspond to mid-1980s fossil fuel combustion since CO
# emissions have held constant (Granier et al., 2011) while CO2 has risen.
#
# All scale factors for years 2000-2013 lumped into a single entry
#==============================================================================
(((CO2CORR

#410 CO2_FOSSFUEL_CORR 0.0489 - - - xy 1 1

## CO2 scale factors from CO oxidation by decade (comment out for now)
##410 CO2_FOSSFUEL_CORR 0.052815/0.054486/0.054902/0.055107/0.053166/0.051611/0.050065/0.048803/0.047060/0.046041 - 1980-1989/1/1/0 C xy 1 1
##410 CO2_FOSSFUEL_CORR 0.045816/0.045153/0.045541/0.045556/0.044800/0.043875/0.042910/0.042206/0.042257/0.042468 - 1990-1999/1/1/0 C xy 1 1
##410 CO2_FOSSFUEL_CORR 0.041495/0.040525/0.040125/0.037853/0.035957/0.034686/0.033538/0.032771/0.031961/0.032118 - 2000-2009/1/1/0 C xy 1 1
##410 CO2_FOSSFUEL_CORR 0.030622/0.029674/0.029029/0.028355                                                       - 2010-2013/1/1/0 C xy 1 1

410 CO2_FOSSFUEL_CORR 0.041495/0.040525/0.040125/0.037853/0.035957/0.034686/0.033538/0.032771/0.031961/0.032118/0.030622/0.029674/0.029029/0.028355 - 2000-2013/1/1/0 C xy 1 1
)))CO2CORR

#==============================================================================
# --- CO2 CONVERSION FACTORS ---
#==============================================================================
420 CH4TOCO2  2.75   - - - xy 1 1
421 CTOCO2    3.6667 - - - xy 1 1
430 NMHCSCALE 0.3333 - - - xy 1 1

#==============================================================================
# ---- TIMES diurnal and weekly scale factors for national fossil fuel CO2 ---
#
# These temporal scale factors are described in Nassar et al. (2013)
#==============================================================================
(((FOSSIL_CDIAC.or.FOSSIL_ODIAC.or.CO2CORR
440 CO2_DIURNAL $ROOT/CO2/v2015-04/FOSSIL/TIMES_diurnal_scale_factors.nc diurnal_scale_factors 2006/1/1/1-24 C xy 1 1
441 CO2_WEEKLY  $ROOT/CO2/v2015-04/FOSSIL/TIMES_weekly_scale_factors.nc  weekly_scale_factors  2006/1/WD/0   C xy 1 1
)))FOSSIL_CDIAC.or.FOSSIL_ODIAC.or.CO2CORR

#==============================================================================
# --- DOMESTIC AVIATION SURFACE CORRECTION FACTOR ---
#
# Regional scale factors slightly less than 1 remove surface contribution from
# domestic aviation in national fossil fuel emisisons, which should instead be
# counted in the air as a 3D field. This scale factor should be applied to the
# main/national fossil fuel emission field if using avaition emissions.
#==============================================================================
(((FOSSIL_CDIAC.or.FOSSIL_ODIAC.or.CO2CORR
480 AVIATION_SURF_CORR $ROOT/CO2/v2015-04/FOSSIL/Aviation_SurfCorr_SclFac.1x1.nc CO2 2004/1/1/0 C xy 1 1
)))FOSSIL_CDIAC.or.FOSSIL_ODIAC.or.CO2CORR

#==============================================================================
# --- CO: HTAP Annual scale factors ---
#==============================================================================
(((HTAP
506 LIQFUEL_THISYR    $ROOT/AnnualScalar/v2014-07/AnnualScalar.geos.1x1.nc COscalar  1985-2010/1/1/0 C xy 1  1
507 LIQFUEL_1985      $ROOT/AnnualScalar/v2014-07/AnnualScalar.geos.1x1.nc COscalar  1985/1/1/0      C xy 1 -1
509 LIQFUEL_2006      $ROOT/AnnualScalar/v2014-07/AnnualScalar.geos.1x1.nc COscalar  2006/1/1/0      C xy 1 -1
510 LIQFUEL_2002      $ROOT/AnnualScalar/v2014-07/AnnualScalar.geos.1x1.nc COscalar  2002/1/1/0      C xy 1 -1
528 LIQFUEL_2008_2010 $ROOT/AnnualScalar/v2014-07/AnnualScalar.geos.1x1.nc COscalar  2008-2010/1/1/0 C xy 1 -1
)))HTAP

)))EMISSIONS

### END SECTION SCALE FACTORS ###

###############################################################################
### BEGIN SECTION MASKS
###############################################################################

# ScalID Name sourceFile sourceVar sourceTime C/R/E SrcDim SrcUnit Oper Lon1/Lat1/Lon2/Lat2

(((EMISSIONS

#==============================================================================
# --- Masks used for CH4 simulation emissions ---
#==============================================================================
(((USE_CH4_DATA

(((Scarpelli_Mexico
1001 MEX_MASK          $ROOT/MASKS/v2018-09/Mexico_Mask.001x001.nc        MASK 2000/1/1/0 C xy 1 1 -118/17/-95/33
1010 MEX_MASK_MIRROR   $ROOT/MASKS/v2018-09/Mexico_Mask_Mirror.001x001.nc MASK 2000/1/1/0 C xy 1 1 -118/17/-95/33
)))Scarpelli_Mexico

(((Scarpelli_Canada
1002 CAN_MASK          $ROOT/MASKS/v2018-09/Canada_Mask.001x001.nc        MASK 2000/1/1/0 C xy 1 1 -141/40/-52/85
1011 CAN_MASK_MIRROR   $ROOT/MASKS/v2018-09/Canada_Mask_Mirror.001x001.nc MASK 2000/1/1/0 C xy 1 1 -141/40/-52/85
)))Scarpelli_Canada

(((GEPA
1008 CONUS_MASK        $ROOT/MASKS/v2018-09/CONUS_Mask.001x001.nc         MASK 2000/1/1/0 C xy 1 1 -140/20/-50/60
1009 CONUS_MASK_MIRROR $ROOT/MASKS/v2018-09/CONUS_Mask_Mirror.001x001.nc  MASK 2000/1/1/0 C xy 1 1 -140/20/-50/60
)))GEPA

)))USE_CH4_DATA

#==============================================================================
# --- Masks used for CO and CO2 simulation emissions ---
#==============================================================================
(((USE_CO_DATA.or.USE_CO2_DATA

(((APEI
2002 CANADA_MASK $ROOT/MASKS/v2018-09/Canada_mask.geos.1x1.nc            MASK     2000/1/1/0 C xy 1 1 -141/40/-52/85
)))APEI

(((NEI2016_MONMEAN
2007 CONUS_MASK  $ROOT/MASKS/v2018-09/CONUS_Mask.01x01.nc                MASK     2000/1/1/0 C xy 1 1 -140/20/-50/60
)))NEI2016_MONMEAN

(((DICE_Africa
2008 AFRICA_MASK $ROOT/MASKS/v2018-09/AF_LANDMASK.geos.05x0666.global.nc LANDMASK 1985/1/1/0 C xy 1 1 -20/-37/54/40
)))DICE_Africa

(((CEDSv2.or.CEDS_GBDMAPS.or.CEDS_GBDMAPS_byFuelType
2009 CHINA_MASK  $ROOT/MASKS/v2018-09/China_mask.generic.1x1.nc          MASK     2000/1/1/0 C xy 1 1  70/10/150/60
)))CEDSv2.or.CEDS_GBDMAPS.or.CEDS_GBDMAPS_byFuelType

)))USE_CO_DATA.or.USE_CO2_DATA

)))EMISSIONS

### END SECTION MASKS ###

###############################################################################
### REFERENCES
###############################################################################
#
# Andres. R.J., J.S Gregg, L. Losey, G. Marland, T.A. Boden (2011), Montly,
# global emissions of carbon dioxide from fossil fuel consumption, Tellus 63B,
# 309-327.
#
# Baker, D.F. et al. (2006) TransCom3 inversion intercomparison: Impact of
# transport model errors on the interannual variability of regional CO2 fluxes,
# 1998-2003, Global Biogeochem., Cy., 20, GB1002, doi:10.29/2004GB002439.
#
# Endresen, O. et al. (2007), A historical reconstuctino of ships fuel
# consuption and emissions, J. Geophys. Res., 112, D12301,
# doi:10.1029/2006JD007630.
#
# Granier, C. et al. (2011), Evolution of anthropogenic and biomass burning
# emissions of air pollutants at global and regional scales during the 1980-2010
# period, Climatic Change, 109:163-190, doi:10.1007/s10584-011-0154-1.
#
# Keller, C.A., M.S. Long, R.M. Yantosca, A.M. DaSilva, S. Pawson, D.J. Jacob
# (2014), HEMCO v1.0: a versatile, ESMF-compliant component for calculating
# emissions in atmospheric models, Geosci., Model Dev., 7, 1409?1417,
# doi:10.5194/gmd-7-1409-2014.
#
# Messerschmidt, J., N. Parazoo, N.M. Deutscher, C. Roehl, T. Warneke, P.O.
# Wennberg, and D. Wunch (2012) Evaluation of atmosphere-biosphere exchange
# estimations with TCCON measurements, Atmos. Chem. Phys. Discussions, 12,
# 12759-12800, doi:10.5194/acpd-12-12759-2012.
#
# Nassar, R., D.B.A. Jones, P. Suntharalingam, J.M. Chen, R.J. Andres, K.J.
# Wecht, R.M. Yantosca, S.S. Kulawik, K.W. Bowman, J.R. Worden, T. Machida and
# H. Matsueda (2010), Modeling global atmospheric CO2 with improved emission
# inventories and CO2 production from the oxidation of other carbon species,
# Geoscientific Model Development, 3, 689-716.
#
# Nassar, R., L. Napier-Linton, K.R. Gurney, R.J. Andres, T. Oda, F.R. Vogel,
# F. Deng (2013), Improving the temporal and spatial distribution of CO2
# emissions from global fossil fuel emission datasets, Journal of Geophysical
# Research: Atmospheres, 118, 917-933, doi:10.1029/2012JD018196.
#
# Oda, T. and S. Maksyutov (2011), A very high-resolution (1 km x 1 km) global
# fossil fuel CO2 emission inventory derived using a point source database and
# satellite observations of nighttime lights, Atmos. Chem. Phys., 11, 543?556,
# doi:10.5194/acp-11-543-2011.
#
# Olsen, S.C. and J.T. Randerson (2004), Differences between surface and column
# atmospheric CO2 and implications for carbon cycle research, J. Geophys. Res.
# 109, D02301, doi:10.1029/2003JD003968.
#
# Olsen, S.C., D.J. Weubbles, B. Owen (2013) Comparison of global 3-D aviation
# datasets, Atmos. Chem. Phys., 13, 429?441, doi:10.5194/acp-13-429-2013.
#
# Simone, N., M. Stettler, S. Eastham, S. Barrett, Aviation Emissions Inventory
# Code (AEIC ) User Manual (R1), Laboratory for Aviation and the Environment,
# Massachusetts Institute of Technology, January 2013, Report No:
# LAE-2013-001-N,
# www.LAE.MIT.edu.
#
# Takahashi, T., et al. (2009), Climatological mean and decadal change in
# surface ocean pCO2, and net sea-air CO2 flux over the global oceans, Deep-Sea
# Res. II, 56(8?10), 554?577, doi:10.1016/j.dsr2.2008.12.009.
#

### END OF HEMCO INPUT FILE ###
#EOC
