Ext_AllowExtrap: .true.
#
PrimaryExports%%
#--------|-------|------|------------|---------------|--------|-------|---------------------|
# Export |       |      |            |___ Refresh ___|____ Factors ___|___ External File ___|
# Name   | Units | Clim |Conservative| Time Template | Offset | Scale | Variable | Template |
#--------|-------|------|------------|---------------|--------|-------|----------|----------|
#
# Notes:
# Units should be in single quotes if they contain whitespace
# Climatology should be Y if the file contains monthly climatology; otherwise it should be N
# Conservative should be Y if units imply mass conservation dependency on regrid method (e.g. value per area)
# Refresh Time Template should be:
#    1. - if the file contains time-invariant constants
#    2. prefixed with F if no time interpolation between data reads (F = fixed)
#    3. 0 if data should be kept up-to-date at all times
#    4. 0:HHMMSS if data occurs in file with frequency HHMMSS and should be updated at that frequency
#    5. %y4-%m2-%d2T%h2:00:00, with any tokens replaced by appropriate constant time value, if multiple time
#       values are in the same file. Note that including an F prefix indicates the data should be treated
#       as constant until the next refresh time. Omission of the F prefix will time-interpolate between
#       the values.
# Minimize whitespace wherever possible to avoid surpassing the ESMF read buffer character limit,
# applicable only for uncommented lines.
#
###############################################################################
###
### Meteorology data
###
###############################################################################
${RUNDIR_MET_EXTDATA_PRIMARY_EXPORTS}
#
#==============================================================================
# --- Fields for lightning emissions (Extension 103) ---
# These fields are stored in State_Met, along with the other met fields
#==============================================================================
FLASH_DENS 1 N Y F0;013000 none none LDENS ./HcoDir/OFFLINE_LIGHTNING/v2020-03/${RUNDIR_MET}/%y4/FLASH_CTH_${RUNDIR_MET}_${RUNDIR_MET_NATIVE_RES}_%y4_%m2.nc4
CONV_DEPTH 1 N Y F0;013000 none none CTH   ./HcoDir/OFFLINE_LIGHTNING/v2020-03/${RUNDIR_MET}/%y4/FLASH_CTH_${RUNDIR_MET}_${RUNDIR_MET_NATIVE_RES}_%y4_%m2.nc4
#
###############################################################################
###
### Land data (not handled by HEMCO)
###
###############################################################################
#
#==============================================================================
# Olson land types
#==============================================================================
# Use conservative fraction regridding to extract land type fraction
OLSON00 1 N F;0  - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON01 1 N F;1  - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON02 1 N F;2  - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON03 1 N F;3  - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON04 1 N F;4  - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON05 1 N F;5  - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON06 1 N F;6  - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON07 1 N F;7  - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON08 1 N F;8  - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON09 1 N F;9  - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON10 1 N F;10 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON11 1 N F;11 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON12 1 N F;12 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON13 1 N F;13 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON14 1 N F;14 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON15 1 N F;15 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON16 1 N F;16 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON17 1 N F;17 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON18 1 N F;18 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON19 1 N F;19 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON20 1 N F;20 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON21 1 N F;21 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON22 1 N F;22 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON23 1 N F;23 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON24 1 N F;24 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON25 1 N F;25 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON26 1 N F;26 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON27 1 N F;27 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON28 1 N F;28 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON29 1 N F;29 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON30 1 N F;30 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON31 1 N F;31 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON32 1 N F;32 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON33 1 N F;33 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON34 1 N F;34 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON35 1 N F;35 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON36 1 N F;36 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON37 1 N F;37 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON38 1 N F;38 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON39 1 N F;39 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON40 1 N F;40 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON41 1 N F;41 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON42 1 N F;42 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON43 1 N F;43 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON44 1 N F;44 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON45 1 N F;45 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON46 1 N F;46 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON47 1 N F;47 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON48 1 N F;48 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON49 1 N F;49 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON50 1 N F;50 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON51 1 N F;51 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON52 1 N F;52 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON53 1 N F;53 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON54 1 N F;54 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON55 1 N F;55 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON56 1 N F;56 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON57 1 N F;57 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON58 1 N F;58 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON59 1 N F;59 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON60 1 N F;60 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON61 1 N F;61 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON62 1 N F;62 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON63 1 N F;63 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON64 1 N F;64 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON65 1 N F;65 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON66 1 N F;66 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON67 1 N F;67 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON68 1 N F;68 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON69 1 N F;69 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON70 1 N F;70 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON71 1 N F;71 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
OLSON72 1 N F;72 - none none OLSON ./ChemDir/Olson_Land_Map_201203/Olson_2001_Land_Map.025x025.generic.nc
#
#==============================================================================
# Yuan-processed MODIS Leaf Area Index
#==============================================================================
#
# Use file with land type stored as 3rd dimension to speed up file read
XLAIMULTI cm2_cm-2 N Y %y4-%m2-%d2T00:00:00 none none XLAIMULTI ./HcoDir/Yuan_XLAI/v2021-06/Condensed_Yuan_proc_MODIS_XLAI.025x025.%y4.nc
#
###############################################################################
###
### HEMCO Base Emissions (Update if HEMCO_Config.rc changes)
###
###############################################################################

#==============================================================================
# %%%%% CH4 emissions %%%%%
#==============================================================================

# --- CH4: Gridded EPA (Maasakkers et al., Environ. Sci. Technol., 2016) ---
GEPA_OIL               molec/cm2/s N    Y -                    none none  emissions_6B_Wastewater_Treatment_Industrial  ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_GAS_PRODUCTION    molec/cm2/s N    Y -                    none none  emissions_1B2b_Natural_Gas_Production         ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_GAS_PROCESSING    molec/cm2/s N    Y -                    none none  emissions_1B2b_Natural_Gas_Processing         ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_GAS_TRANSMISSION  molec/cm2/s N    Y -                    none none  emissions_1B2b_Natural_Gas_Transmission       ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_GAS_DISTRIBUTION  molec/cm2/s N    Y -                    none none  emissions_1B2b_Natural_Gas_Distribution       ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_COAL_UNDERGROUND  molec/cm2/s N    Y -                    none none  emissions_1B1a_Coal_Mining_Underground        ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_COAL_SURFACE      molec/cm2/s N    Y -                    none none  emissions_1B1a_Coal_Mining_Surface            ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_COAL_ABANDONED    molec/cm2/s N    Y -                    none none  emissions_1B1a_Abandoned_Coal                 ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_LIVESTOCK__4A     molec/cm2/s N    Y -                    none none  emissions_4A_Enteric_Fermentation             ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_LIVESTOCK__4B     molec/cm2/s N    Y -                    none none  emissions_4B_Manure_Management                ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_LANDFILLS_MUNI    molec/cm2/s N    Y -                    none none  emissions_6A_Landfills_Municipal              ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_LANDFILLS_INDU    molec/cm2/s N    Y -                    none none  emissions_6A_Landfills_Industrial             ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_WASTEWATER_DOME   molec/cm2/s N    Y -                    none none  emissions_6B_Wastewater_Treatment_Domestic    ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_WASTEWATER_INDU   molec/cm2/s N    Y -                    none none  emissions_6B_Wastewater_Treatment_Industrial  ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_RICE              molec/cm2/s N    Y -                    none none  emissions_2C2_Ferroalloy_Production           ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_OTHER__1A_M       molec/cm2/s N    Y -                    none none  emissions_1A_Combustion_Mobile                ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_OTHER__1A_S       molec/cm2/s N    Y -                    none none  emissions_1A_Combustion_Stationary            ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_OTHER__2B5        molec/cm2/s N    Y -                    none none  emissions_2B5_Petrochemical_Production        ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_OTHER__2C2        molec/cm2/s N    Y -                    none none  emissions_2C2_Ferroalloy_Production           ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_OTHER__4F         molec/cm2/s 2012 Y F%y4-%m2-01T00:00:00 none none  emissions_4F_Field_Burning                    ./HcoDir/CH4/v2022-11/GEPA/GEPA_Monthly.nc
#GEPA_OTHER__5          molec/cm2/s 2012 Y F%y4-%m2-%d2T00:00:00 none none emissions_5_Forest_Fires                      ./HcoDir/CH4/v2022-11/GEPA/GEPA_Daily.nc
GEPA_OTHER__6D         molec/cm2/s N    Y -                    none none  emissions_6D_Composting                       ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_COAST_OIL         molec/cm2/s N    Y -                    none none  emissions_1B2a_Petroleum                      ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_COAST_GAS_PR      molec/cm2/s N    Y -                    none none  emissions_1B2b_Natural_Gas_Production         ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_COAST_GAS_PC      molec/cm2/s N    Y -                    none none  emissions_1B2b_Natural_Gas_Processing         ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_COAST_GAS_TR      molec/cm2/s N    Y -                    none none  emissions_1B2b_Natural_Gas_Transmission       ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_COAST_GAS_DS      molec/cm2/s N    Y -                    none none  emissions_1B2b_Natural_Gas_Distribution       ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_COAST_COAL_U      molec/cm2/s N    Y -                    none none  emissions_1B1a_Coal_Mining_Underground        ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_COAST_COAL_S      molec/cm2/s N    Y -                    none none  emissions_1B1a_Coal_Mining_Surface            ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_COAST_COAL_A      molec/cm2/s N    Y -                    none none  emissions_1B1a_Abandoned_Coal                 ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_COAST_LVSTK_F     molec/cm2/s N    Y -                    none none  emissions_4A_Enteric_Fermentation             ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_COAST_LVSTK_M     molec/cm2/s N    Y -                    none none  emissions_4B_Manure_Management                ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_COAST_LDF_M       molec/cm2/s N    Y -                    none none  emissions_6A_Landfills_Municipal              ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_COAST_LDF_I       molec/cm2/s N    Y -                    none none  emissions_6A_Landfills_Industrial             ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_COAST_WSTW_D      molec/cm2/s N    Y -                    none none  emissions_6B_Wastewater_Treatment_Domestic    ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_COAST_WSTW_I      molec/cm2/s N    Y -                    none none  emissions_6B_Wastewater_Treatment_Industrial  ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_COAST_RICE        molec/cm2/s N    Y -                    none none  emissions_4C_Rice_Cultivation                 ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_COAST_OTH_1AM     molec/cm2/s N    Y -                    none none  emissions_1A_Combustion_Mobile                ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_COAST_OTH_1AS     molec/cm2/s N    Y -                    none none  emissions_1A_Combustion_Stationary            ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_COAST_OTH_2B5     molec/cm2/s N    Y -                    none none  emissions_2B5_Petrochemical_Production        ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_COAST_OTH_2C2     molec/cm2/s N    Y -                    none none  emissions_2C2_Ferroalloy_Production           ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc
GEPA_COAST_OTH_4F      molec/cm2/s 2012 Y F%y4-%m2-01T00:00:00 none none  emissions_4F_Field_Burning                    ./HcoDir/CH4/v2022-11/GEPA/GEPA_Monthly.nc
#GEPA_COAST_OTH_5       molec/cm2/s 2012 Y F%y4-%m2-%d2T00:00:00 none none emissions_5_Forest_Fires                       ./HcoDir/CH4/v2022-11/GEPA/GEPA_Daily.nc
GEPA_COAST_OTH_6D      molec/cm2/s N    Y -                    none none  emissions_6D_Composting                       ./HcoDir/CH4/v2022-11/GEPA/GEPA_Annual.nc

# --- Scarpelli Mexico ---
MEX_OIL                      kg/m2/s N Y - none none  emis_ch4  ./HcoDir/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_oil_2015.nc
MEX_GAS                      kg/m2/s N Y - none none  emis_ch4  ./HcoDir/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_gas_2015.nc
MEX_COAL                     kg/m2/s N Y - none none  emis_ch4  ./HcoDir/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_coal_2015.nc
MEX_LIVESTOCK_A              kg/m2/s N Y - none none  emis_ch4  ./HcoDir/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_livestock_A_2015.nc
MEX_LIVESTOCK_B              kg/m2/s N Y - none none  emis_ch4  ./HcoDir/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_livestock_B_2015.nc
MEX_LANDFILLS                kg/m2/s N Y - none none  emis_ch4  ./HcoDir/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_landfill_2015.nc
MEX_WASTEWATER               kg/m2/s N Y - none none  emis_ch4  ./HcoDir/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_waste_2015.nc
MEX_RICE                     kg/m2/s N Y - none none  emis_ch4  ./HcoDir/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_rice_2015.nc
MEX_OTHER                    kg/m2/s N Y - none none  emis_ch4  ./HcoDir/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_other_anthro_2015.nc
MEX_OIL_COAST                kg/m2/s N Y - none none  emis_ch4  ./HcoDir/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_oil_2015.nc
MEX_GAS_COAST                kg/m2/s N Y - none none  emis_ch4  ./HcoDir/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_gas_2015.nc
MEX_COAL_COAST               kg/m2/s N Y - none none  emis_ch4  ./HcoDir/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_coal_2015.nc
MEX_LIVESTOCK_A_COAST        kg/m2/s N Y - none none  emis_ch4  ./HcoDir/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_livestock_A_2015.nc
MEX_LIVESTOCK_B_COAST        kg/m2/s N Y - none none  emis_ch4  ./HcoDir/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_livestock_B_2015.nc
MEX_LANDFILLS_COAST          kg/m2/s N Y - none none  emis_ch4  ./HcoDir/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_landfill_2015.nc
MEX_WASTEWATER_COAST         kg/m2/s N Y - none none  emis_ch4  ./HcoDir/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_waste_2015.nc
MEX_RICE_COAST               kg/m2/s N Y - none none  emis_ch4  ./HcoDir/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_rice_2015.nc
MEX_OTHER_COAST              kg/m2/s N Y - none none  emis_ch4  ./HcoDir/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_other_anthro_2015.nc

# --- Scarpelli Canada ---
CAN_OIL_GAS_COMBUSTION       kg/m2/s N Y - none none  oil_gas_combustion_total   ./HcoDir/CH4/v2022-01/Scarpelli_Canada/can_emis_oil_gas_combustion_2018.nc
CAN_OIL_GAS_LEAKAGE          kg/m2/s N Y - none none  oil_gas_leakage_total      ./HcoDir/CH4/v2022-01/Scarpelli_Canada/can_emis_oil_gas_leakage_2018.nc
CAN_OIL_GAS_VENT_FLARE       kg/m2/s N Y - none none  oil_gas_vent_flare_total   ./HcoDir/CH4/v2022-01/Scarpelli_Canada/can_emis_oil_gas_vent_flare_2018.nc
CAN_COAL                     kg/m2/s N Y - none none  coal_total                 ./HcoDir/CH4/v2022-01/Scarpelli_Canada/can_emis_coal_2018.nc
CAN_LIVESTOCK                kg/m2/s N Y - none none  livestock_total            ./HcoDir/CH4/v2022-01/Scarpelli_Canada/can_emis_livestock_2018.nc
CAN_SOLID_WASTE              kg/m2/s N Y - none none  solid_waste_total          ./HcoDir/CH4/v2022-01/Scarpelli_Canada/can_emis_solid_waste_2018.nc
CAN_WASTEWATER               kg/m2/s N Y - none none  wastewater_total           ./HcoDir/CH4/v2022-01/Scarpelli_Canada/can_emis_wastewater_2018.nc
CAN_OTHER                    kg/m2/s N Y - none none  other_minor_sources_total  ./HcoDir/CH4/v2022-01/Scarpelli_Canada/can_emis_other_minor_sources_2018.nc
CAN_OIL_GAS_COMBUSTION_COAST kg/m2/s N Y - none none  oil_gas_combustion_total   ./HcoDir/CH4/v2022-01/Scarpelli_Canada/can_emis_oil_gas_combustion_2018.nc
CAN_OIL_GAS_LEAKAGE_COAST    kg/m2/s N Y - none none  oil_gas_leakage_total      ./HcoDir/CH4/v2022-01/Scarpelli_Canada/can_emis_oil_gas_leakage_2018.nc
CAN_OIL_GAS_VENT_FLARE_COAST kg/m2/s N Y - none none  oil_gas_vent_flare_total   ./HcoDir/CH4/v2022-01/Scarpelli_Canada/can_emis_oil_gas_vent_flare_2018.nc
CAN_COAL_COAST               kg/m2/s N Y - none none  coal_total                 ./HcoDir/CH4/v2022-01/Scarpelli_Canada/can_emis_coal_2018.nc
CAN_LIVESTOCK_COAST          kg/m2/s N Y - none none  livestock_total            ./HcoDir/CH4/v2022-01/Scarpelli_Canada/can_emis_livestock_2018.nc
CAN_SOLID_WASTE_COAST        kg/m2/s N Y - none none  solid_waste_total          ./HcoDir/CH4/v2022-01/Scarpelli_Canada/can_emis_solid_waste_2018.nc
CAN_WASTEWATER_COAST         kg/m2/s N Y - none none  wastewater_total           ./HcoDir/CH4/v2022-01/Scarpelli_Canada/can_emis_wastewater_2018.nc
CAN_OTHER_COAST              kg/m2/s N Y - none none  other_minor_sources_total  ./HcoDir/CH4/v2022-01/Scarpelli_Canada/can_emis_other_minor_sources_2018.nc

# --- GFEI ---
GFEI_CH4_OIL  kg/m2/s N Y - none none  emis_ch4  ./HcoDir/CH4/v2022-01/GFEIv2/Global_Fuel_Exploitation_Inventory_v2_2019_Oil_All.nc
GFEI_CH4_GAS  kg/m2/s N Y - none none  emis_ch4  ./HcoDir/CH4/v2022-01/GFEIv2/Global_Fuel_Exploitation_Inventory_v2_2019_Gas_All.nc
GFEI_CH4_COAL kg/m2/s N Y - none none  emis_ch4  ./HcoDir/CH4/v2022-01/GFEIv2/Global_Fuel_Exploitation_Inventory_v2_2019_Coal.nc

# --- EDGAR v6.0 ---
CH4_OIL__1B2a          kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2022-11/EDGARv6/%y4/v6.0_CH4_%y4_PRO_OIL.0.1x0.1.nc
CH4_OIL__1B2c          kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2022-11/EDGARv6/%y4/v6.0_CH4_%y4_PRO_GAS.0.1x0.1.nc
CH4_COAL__1B1a         kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2022-11/EDGARv6/%y4/v6.0_CH4_%y4_PRO_COAL.0.1x0.1.nc
CH4_LIVESTOCK__4A      kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2022-11/EDGARv6/%y4/v6.0_CH4_%y4_ENF.0.1x0.1.nc
CH4_LIVESTOCK__4B      kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2022-11/EDGARv6/%y4/v6.0_CH4_%y4_MNM.0.1x0.1.nc
CH4_LANDFILLS__6A_6D   kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2022-11/EDGARv6/%y4/v6.0_CH4_%y4_SWD_LDF.0.1x0.1.nc
CH4_WASTEWATER__6B     kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2022-11/EDGARv6/%y4/v6.0_CH4_%y4_WWT.0.1x0.1.nc
CH4_RICE__4C_4D        kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2022-11/EDGARv6/%y4/v6.0_CH4_%y4_AGS.0.1x0.1.nc
CH4_OTHER__1A1_1B1_1B2 kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2022-11/EDGARv6/%y4/v6.0_CH4_%y4_REF_TRF.0.1x0.1.nc
CH4_OTHER__1A1a        kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2022-11/EDGARv6/%y4/v6.0_CH4_%y4_ENE.0.1x0.1.nc
CH4_OTHER__1A2         kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2022-11/EDGARv6/%y4/v6.0_CH4_%y4_IND.0.1x0.1.nc
CH4_OTHER__1A3a_CDS    kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2022-11/EDGARv6/%y4/v6.0_CH4_%y4_TNR_Aviation_CDS.0.1x0.1.nc
CH4_OTHER__1A3a_CRS    kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2022-11/EDGARv6/%y4/v6.0_CH4_%y4_TNR_Aviation_CRS.0.1x0.1.nc
CH4_OTHER__1A3a_LTO    kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2022-11/EDGARv6/%y4/v6.0_CH4_%y4_TNR_Aviation_LTO.0.1x0.1.nc
CH4_OTHER__1A3b        kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2022-11/EDGARv6/%y4/v6.0_CH4_%y4_TRO_noRES.0.1x0.1.nc
CH4_OTHER__1A3c_1A3e   kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2022-11/EDGARv6/%y4/v6.0_CH4_%y4_TNR_Other.0.1x0.1.nc
CH4_OTHER__1A3d_1C2    kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2022-11/EDGARv6/%y4/v6.0_CH4_%y4_TNR_Ship.0.1x0.1.nc
CH4_OTHER__1A4         kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2022-11/EDGARv6/%y4/v6.0_CH4_%y4_RCO.0.1x0.1.nc
CH4_OTHER__2B          kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2022-11/EDGARv6/%y4/v6.0_CH4_%y4_CHE.0.1x0.1.nc
CH4_OTHER__2C          kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2022-11/EDGARv6/%y4/v6.0_CH4_%y4_IRO.0.1x0.1.nc
CH4_OTHER__4F          kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2022-11/EDGARv6/%y4/v6.0_CH4_%y4_AWB.0.1x0.1.nc
CH4_OTHER__6C          kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2022-11/EDGARv6/%y4/v6.0_CH4_%y4_SWD_INC.0.1x0.1.nc

# --- CH4: CEDS (historical) or Shared Socioeconomic Pathways (future) ---
# --- NOTE: This is only for GCAP2 meteorology, so we can comment out  ---
##CMIP6_CH4_AGR  kg/m2/s N Y %y4-%m2-01T00:00:00  none none  CH4_agr  ./HcoDir/CMIP6/v2021-01/$GCAP2SCENARIO/$GCAP2SCENARIO_%y4.nc4
##CMIP6_CH4_ENE  kg/m2/s N Y %y4-%m2-01T00:00:00  none none  CH4_ene  ./HcoDir/CMIP6/v2021-01/$GCAP2SCENARIO/$GCAP2SCENARIO_%y4.nc4
##CMIP6_CH4_IND  kg/m2/s N Y %y4-%m2-01T00:00:00  none none  CH4_ind  ./HcoDir/CMIP6/v2021-01/$GCAP2SCENARIO/$GCAP2SCENARIO_%y4.nc4
##CMIP6_CH4_TRA  kg/m2/s N Y %y4-%m2-01T00:00:00  none none  CH4_tra  ./HcoDir/CMIP6/v2021-01/$GCAP2SCENARIO/$GCAP2SCENARIO_%y4.nc4
##CMIP6_CH4_RCO  kg/m2/s N Y %y4-%m2-01T00:00:00  none none  CH4_rco  ./HcoDir/CMIP6/v2021-01/$GCAP2SCENARIO/$GCAP2SCENARIO_%y4.nc4
##CMIP6_CH4_SLV  kg/m2/s N Y %y4-%m2-01T00:00:00  none none  CH4_slv  ./HcoDir/CMIP6/v2021-01/$GCAP2SCENARIO/$GCAP2SCENARIO_%y4.nc4
##CMIP6_CH4_WST  kg/m2/s N Y %y4-%m2-01T00:00:00  none none  CH4_wst  ./HcoDir/CMIP6/v2021-01/$GCAP2SCENARIO/$GCAP2SCENARIO_%y4.nc4
##CMIP6_CH4_SHP  kg/m2/s N Y %y4-%m2-01T00:00:00  none none  CH4_shp  ./HcoDir/CMIP6/v2021-01/$GCAP2SCENARIO/$GCAP2SCENARIO_%y4.nc4

# --- CH4: BB4MIPs historical / SSP future biomass burning inventories ---
# CMIP6_BB_CH4 %y4-%m2-01T00:00:00 none none  CH4_shp  ./HcoDir/CMIP6/v2021-01/$GCAP2SCENARIO/$GCAP2SCENARIO_%y4.nc4

# --- CH4: Updated GFED4 biomass burning (Yuzhong Zhang) ---
UPDATED_GFED4_CH4 kg/m2/s N Y F%y4-%m2-%d2T00:00:00 none none  CH4 ./HcoDir/CH4/v2020-09/GFED/GFED4_CH4_FIRE_%y4.nc

# --- JPL WetCHARTs v1.0                                    ---
# --- (Bloom et al., https://doi.org/10.3334/ORNLDAAC/1502) ---
JPLW_CH4  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none  emi_ch4  ./HcoDir/CH4/v2020-09/JPL_WetCharts/HEensemble/JPL_WetCharts_2010-2019.Ensemble_Mean.0.5x0.5.nc

# --- Geological seeps ---
CH4_SEEPS kg/m2/s N Y - none none  emi_ch4  ./HcoDir/CH4/v2020-04/Seeps/Etiope_CH4GeologicalEmis_ScaledToHmiel.1x1.nc

# --- Emissions from lakes (Maasakkers et al.) ---
# --- NOTE: Superseded by WetCHARTs            ---
#CH4_LAKES kg/m2/s N Y F%y4-%m2-01T00:00:00 none none  emi_ch4  ./HcoDir/CH4/v2022-11/Lakes/Maasakkers_Lakes_%y4.01x01.nc

# --- Emissions from termites (Fung et al 1991) ---
CH4_TERMITES    kg/m2/s N Y - none none  CH4  ./HcoDir/CH4/v2022-11/4x5/termites.geos.4x5.nc

# --- Soil absoprtion (Fung et al 1991) ---
# --- NOTE: Superseded by MeMo          ---
#CH4_SOILABSORB  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none  CH4  ./HcoDir/CH4/v2019-10/Fung_SoilAbs/Soil_Absorption_4x5_%y4.nc

# --- Soil absorption fro MeMo model (Murguia-Flores et al. 2018, GMD) ---
# --- Note; Use climatology by default                                 ---
CH4_SOILABSORB kg/m2/s  2009 Y F%y4-%m2-01T00:00:00 none none  CH4uptake  ./HcoDir/CH4/v2019-10/MeMo_SoilAbs/MeMo_CH4uptake_Climatology.nc
#CH4_SOILABSORB kg/m2/s N    Y F%y4-%m2-01T00:00:00 none none  CH4uptake  ./HcoDir/CH4/v2019-10/MeMo_SoilAbs/MeMo_CH4uptake_1990-2009.nc

# --- RCP future emission scenarios ---
#RCP3PD_CH4 kg/m2/s N Y %y4-01-01T00:00:00T00:00:00 none none CO ACCMIP  ./HcoDir/v2020-07/RCP_3PD/RCPs_anthro_CH4_2005-2100_23474.nc
#RCP45_CH4  kg/m2/s N Y %y4-01-01T00:00:00T00:00:00 none none CO ACCMIP  ./HcoDir/v2020-07/RCP_45/RCPs_anthro_CH4_2005-2100_27424.nc
#RCP60_CH4  kg/m2/s N Y %y4-01-01T00:00:00T00:00:00 none none CO ACCMIP  ./HcoDir/v2020-07/RCP_60/RCPs_anthro_CH4_2005-2100_43190.nc
#RCP85_CH4  kg/m2/s N Y %y4-01-01T00:00:00T00:00:00 none none CO ACCMIP  ./HcoDir/v2020-07/RCP_85/RCPs_anthro_CH4_2005-2100_43533.nc

#==============================================================================
# %%%%% CO emissions %%%%%
#==============================================================================

# --- CEDSv2 emissions         ---
# --- Default global inventory ---
CEDS_CO_AGR  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO_agr  ./HcoDir/CEDS/v2021-06/%y4/CO-em-anthro_CMIP_CEDS_%y4.nc
CEDS_CO_ENE  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO_ene  ./HcoDir/CEDS/v2021-06/%y4/CO-em-anthro_CMIP_CEDS_%y4.nc
CEDS_CO_IND  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO_ind  ./HcoDir/CEDS/v2021-06/%y4/CO-em-anthro_CMIP_CEDS_%y4.nc
CEDS_CO_TRA  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO_tra  ./HcoDir/CEDS/v2021-06/%y4/CO-em-anthro_CMIP_CEDS_%y4.nc
CEDS_CO_RCO  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO_rco  ./HcoDir/CEDS/v2021-06/%y4/CO-em-anthro_CMIP_CEDS_%y4.nc
CEDS_CO_SLV  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO_slv  ./HcoDir/CEDS/v2021-06/%y4/CO-em-anthro_CMIP_CEDS_%y4.nc
CEDS_CO_WST  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO_wst  ./HcoDir/CEDS/v2021-06/%y4/CO-em-anthro_CMIP_CEDS_%y4.nc

# --- CEDSv2 ship emisisons ---
CEDS_CO_SHP  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO_shp  ./HcoDir/CEDS/v2021-06/%y4/CO-em-anthro_CMIP_CEDS_%y4.nc

# --- APEI Canada (APEI) ---
#APEI_CO kg/m2/s N Y F%y4-01-01T00:00:00 none none CO  ./HcoDir/APEI/v2016-11/APEI.0.1x0.1.nc

# --- EPA NEI2016 USA monthly mean (NEI2016_MONMEAN) ---
#EPA16_CO__airportsCO    kg/m2/s  Y Y F2016-%m2-01T00:00:00 none none CO  ./HcoDir/NEI2016/v2021-06/2016fh_16j_airports_0pt1degree_month_%m2.ncf
#EPA16_CO__nonptCO       kg/m2/s  Y Y F2016-%m2-01T00:00:00 none none CO  ./HcoDir/NEI2016/v2021-06/2016fh_16j_nonpt_0pt1degree_month_%m2.ncf
#EPA16_CO__nonroadCO     kg/m2/s  Y Y F2016-%m2-01T00:00:00 none none CO  ./HcoDir/NEI2016/v2021-06/2016fh_16j_nonroad_0pt1degree_month_%m2.ncf
#EPA16_CO__npogCO        kg/m2/s  Y Y F2016-%m2-01T00:00:00 none none CO  ./HcoDir/NEI2016/v2021-06/2016fh_16j_np_oilgas_0pt1degree_month_%m2.ncf
#EPA16_CO__onroadCO      kg/m2/s  Y Y F2016-%m2-01T00:00:00 none none CO  ./HcoDir/NEI2016/v2021-06/2016fh_16j_onroad_0pt1degree_month_%m2.ncf
#EPA16_CO__onroad_caCO   kg/m2/s  Y Y F2016-%m2-01T00:00:00 none none CO  ./HcoDir/NEI2016/v2021-06/2016fh_16j_onroad_ca_adj_0pt1degree_month_%m2.ncf
#EPA16_CO__railCO        kg/m2/s  Y Y F2016-%m2-01T00:00:00 none none CO  ./HcoDir/NEI2016/v2021-06/2016fh_16j_rail_0pt1degree_month_%m2.ncf
#EPA16_CO__rwcCO         kg/m2/s  Y Y F2016-%m2-01T00:00:00 none none CO  ./HcoDir/NEI2016/v2021-06/2016fh_16j_rwc_0pt1degree_month_%m2.ncf
#EPA16_CO__c1c2CO        kg/m2/s  Y Y F2016-%m2-01T00:00:00 none none CO  ./HcoDir/NEI2016/v2021-06/2016fh_16j_emln_cmv_c1c2_12_0pt1degree_3D_month_%m2.ncf
#EPA16_CO__c3CO          kg/m2/s  Y Y F2016-%m2-01T00:00:00 none none CO  ./HcoDir/NEI2016/v2021-06/2016fh_16j_emln_cmv_c3_12_0pt1degree_3D_month_%m2.ncf
#EPA16_CO__pteguCO       kg/m2/s  Y Y F2016-%m2-01T00:00:00 none none CO  ./HcoDir/NEI2016/v2021-06/2016fh_16j_emln_ptegu_0pt1degree_3D_month_%m2.ncf
#EPA16_CO__ptogCO        kg/m2/s  Y Y F2016-%m2-01T00:00:00 none none CO  ./HcoDir/NEI2016/v2021-06/2016fh_16j_emln_pt_oilgas_allinln_0pt1degree_3D_month_%m2.ncf
#EPA16_CO__ptnonipmCO    kg/m2/s  Y Y F2016-%m2-01T00:00:00 none none CO  ./HcoDir/NEI2016/v2021-06/2016fh_16j_emln_ptnonipm_allinln_0pt1degree_3D_month_%m2.ncf
#EPA16_CO__onroad_canCO  kg/m2/s  Y Y F2016-%m2-01T00:00:00 none none CO  ./HcoDir/NEI2016/v2021-06/2016fh_16j_onroad_can_0pt1degree_month_%m2.ncf
#EPA16_CO__onroad_mexCO  kg/m2/s  Y Y F2016-%m2-01T00:00:00 none none CO  ./HcoDir/NEI2016/v2021-06/2016fh_16j_onroad_mex_0pt1degree_month_%m2.ncf
#EPA16_CO__otharCO       kg/m2/s  Y Y F2016-%m2-01T00:00:00 none none CO  ./HcoDir/NEI2016/v2021-06/2016fh_16j_othar_0pt1degree_month_%m2.ncf
#EPA16_CO__othptCO       kg/m2/s  Y Y F2016-%m2-01T00:00:00 none none CO  ./HcoDir/NEI2016/v2021-06/2016fh_16j_emln_othpt_0pt1degree_month_%m2.ncf

# --- DICE-Africa emission inventory (DICE_Africa) ---
## Scale factors convert from g/m2/yr to kg/m2/s
#DICE_CARS_CO          g/m2/yr N Y - none 3.171e-11 CO  ./HcoDir/DICE_Africa/v2016-10/DICE-Africa-cars-2013-v01-4Oct2016.nc
#DICE_MOTORCYCLES_CO   g/m2/yr N Y - none 3.171e-11 CO  ./HcoDir/DICE_Africa/v2016-10/DICE-Africa-motorcycles-2013-v01-4Oct2016.nc
#DICE_BACKUPGEN_CO     g/m2/yr N Y - none 3.171e-11 CO  ./HcoDir/DICE_Africa/v2016-10/DICE-Africa-generator-use-2013-v01-4Oct2016.nc
#DICE_CHARCOALPROD_CO  g/m2/yr N Y - none 3.171e-11 CO  ./HcoDir/DICE_Africa/v2016-10/DICE-Africa-charcoal-production-2013-v01-4Oct2016.nc
#DICE_GASFLARE_CO      g/m2/yr N Y - none 3.171e-11 CO  ./HcoDir/DICE_Africa/v2016-10/DICE-Africa-gas-flares-2013-v01-4Oct2016.nc
#DICE_AGBURNING_CO     g/m2/yr N Y - none 3.171e-11 CO  ./HcoDir/DICE_Africa/v2016-10/DICE-Africa-household-crop-residue-use-2013-v01-4Oct2016.nc
#DICE_CHARCOALUSE_CO   g/m2/yr N Y - none 3.171e-11 CO  ./HcoDir/DICE_Africa/v2016-10/DICE-Africa-charcoal-use-2013-v01-4Oct2016.nc
#DICE_KEROSENE_CO      g/m2/yr N Y - none 3.171e-11 CO  ./HcoDir/DICE_Africa/v2016-10/DICE-Africa-kerosene-use-2013-v01-4Oct2016.nc
#DICE_OILREFINING_CO   g/m2/yr N Y - none 3.171e-11 CO  ./HcoDir/DICE_Africa/v2016-10/DICE-Africa-adhoc-oil-refining-2006-v01-4Oct2016.nc
#DICE_HOUSEFUELWOOD_CO g/m2/yr N Y - none 3.171e-11 CO  ./HcoDir/DICE_Africa/v2016-10/DICE-Africa-household-fuelwood-use-2013-v01-4Oct2016.nc
#DICE_OTHERFUELWOOD_CO g/m2/yr N Y - none 3.171e-11 CO  ./HcoDir/DICE_Africa/v2016-10/DICE-Africa-other-fuelwood-use-2013-v01-4Oct2016.nc
##
##  --- Efficient Combusions Emissions from EDGARv43    ---
##  --- This covers sectors not included in DICE-Africa ---
#AF_EDGAR_CO_POW   kg/m2/s N Y F%y4-01-01T00:00:00 none none emi_co  ./HcoDir/EDGARv43/v2016-11/EDGAR_v43.CO.POW.0.1x0.1.nc
#AF_EDGAR_CO_ENG   kg/m2/s N Y F%y4-01-01T00:00:00 none none emi_co  ./HcoDir/EDGARv43/v2016-11/EDGAR_v43.CO.ENG.0.1x0.1.nc
#AF_EDGAR_CO_IND   kg/m2/s N Y F%y4-01-01T00:00:00 none none emi_co  ./HcoDir/EDGARv43/v2016-11/EDGAR_v43.CO.IND.0.1x0.1.nc
#AF_EDGAR_CO_TNG   kg/m2/s N Y F%y4-01-01T00:00:00 none none emi_co  ./HcoDir/EDGARv43/v2016-11/EDGAR_v43.CO.TNG.0.1x0.1.nc
#AF_EDGAR_CO_PPA   kg/m2/s N Y F%y4-01-01T00:00:00 none none emi_co  ./HcoDir/EDGARv43/v2016-11/EDGAR_v43.CO.PPA.0.1x0.1.nc
#AF_EDGAR_CO_SWD   kg/m2/s N Y F%y4-01-01T00:00:00 none none emi_co  ./HcoDir/EDGARv43/v2016-11/EDGAR_v43.CO.SWD.0.1x0.1.nc

# --- HTAP v2 (HTAP) ---
#HTAP_CO_IND    kg/m2/s N Y F%y4-01-01T00:00:00 none none emi_co  ./HcoDir/HTAP/v2015-03/CO/EDGAR_HTAP_CO_INDUSTRY.generic.01x01.nc
#HTAP_CO_POW    kg/m2/s N Y F%y4-01-01T00:00:00 none none emi_co  ./HcoDir/HTAP/v2015-03/CO/EDGAR_HTAP_CO_ENERGY.generic.01x01.nc
#HTAP_CO_RES    kg/m2/s N Y F%y4-01-01T00:00:00 none none emi_co  ./HcoDir/HTAP/v2015-03/CO/EDGAR_HTAP_CO_RESIDENTIAL.generic.01x01.nc
#HTAP_CO_TRA    kg/m2/s N Y F%y4-01-01T00:00:00 none none emi_co  ./HcoDir/HTAP/v2015-03/CO/EDGAR_HTAP_CO_TRANSPORT.generic.01x01.nc
##HTAP_CO_AIR1  kg/m2/s N Y F%y4-01-01T00:00:00 none none emi_co  ./HcoDir/HTAP/v2015-03/CO/EDGAR_HTAP_CO_AIR_LTO.generic.01x01.nc
##HTAP_CO_AIR2  kg/m2/s N Y F%y4-01-01T00:00:00 none none emi_co  ./HcoDir/HTAP/v2015-03/CO/EDGAR_HTAP_CO_AIR_CDS.generic.01x01.nc
##HTAP_CO_AIR3  kg/m2/s N Y F%y4-01-01T00:00:00 none none emi_co  ./HcoDir/HTAP/v2015-03/CO/EDGAR_HTAP_CO_AIR_CRS.generic.01x01.nc

# --- AEIC 2019 aircraft (AEIC) ---
#AEIC19_DAILY_CO   kg/m2/s 2019 Y F%y4-%m2-%d2T00:00:00 none none CO ./HcoDir/AEIC2019/v2022-03/2019/%m2/AEIC_2019%m2%d2.0.5x0.625.36L.nc
AEIC19_MONMEAN_CO kg/m2/s 2019 Y F%y4-%m2-%d2T00:00:00 none none CO ./HcoDir/AEIC2019/v2022-03/2019_monmean/AEIC_monmean_2019%m2.0.5x0.625.36L.nc

# --- RCP future emissions scenarios ---
#RCP3PD_CO   kg/m2/s N Y %y4-01-01T00:00:00 none none ACCMIP ./HcoDir/RCP/v2020-07/RCP_3PD/RCPs_anthro_CO_2005-2100_23474.nc
#RCP45_CO   kg/m2/s N Y %y4-01-01T00:00:00  none none ACCMIP ./HcoDir/RCP/v2020-07/RCP_45/RCPs_anthro_CO_2005-2100_27424.nc
#RCP60_CO   kg/m2/s N Y %y4-01-01T00:00:00  none none ACCMIP ./HcoDir/RCP/v2020-07/RCP_60/RCPs_anthro_CO_2005-2100_43190.nc
#RCP85_CO   kg/m2/s N Y %y4-01-01T00:00:00  none none ACCMIP ./HcoDir/RCP/v2020-07/RCP_85/RCPs_anthro_CO_2005-2100_43533.nc

#==============================================================================
# --- CO2 emissions ---
#==============================================================================

# --- CO2 Fossil Fuel ---
#FOSSILCO2_CDIAC  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO2 ./HcoDir/CO2/v2014-09/FOSSIL/CDIAC_v2014.monthly.generic.1x1.nc
FOSSILCO2_ODIAC    kg/m2/s N Y  F%y4-%m2-01T00:00:00 none none CO2 ./HcoDir/CO2/v2022-11/FOSSIL/ODIAC_CO2.monthly.generic.1x1.nc

# --- CO2 Ocean Exchange ---
#OCEANCO2_TAKA_ANNUAL   kg/m2/s N Y - none none CO2 ./HcoDir/CO2/v2022-11/OCEAN/Taka2009_CO2_Annual.nc
OCEANCO2_TAKA_MONTHLY  kg/m2/s 2000 Y F%y4-%m2-01T00:00:00 none none CO2 ./HcoDir/CO2/v2022-11/OCEAN/Taka2009_CO2_Monthly.nc

# --- Balanced biosphere exchange ---
#BBIOCO2_DIURNAL kg/m2/s N Y F%y4-%m2-%d2T%h2-01-01 none none CO2 ./HcoDir/CO2/v2014-09/BBIO/BBIO_diurnal_CO2.nc
SIB_BBIO_CO2 kg/m2/s N Y F%y4-%m2-%d2T%h2-01-01 none none CO2 ./HcoDir/CO2/v2022-11/BIO/SiB3_3hr_NEP.nc

# --- Net Terrestrial Exchange ---
CO2_NET_TERRESTRIAL kg/m2/s N Y - none none CO2 ./HcoDir/CO2/v2022-11/BIO/Net_terrestrial_exch_5.29Pg.generic.1x1.nc

# --- CO2 Ship Emissions ---
#ICOADS_CO2_SHIP kg/m2/s N Y F2004-%m2-01T00:00:00 none none CO2 ./HcoDir/ICOADS_SHIP/v2014-07/ICOADS_ship_CO2_2004.generic.1x1.nc
CEDS_CO2_SHP kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO2_shp ./HcoDir/CEDS/v2018-08/%y4/CO2-em-anthro_CMIP_CEDS_%y4.nc

# --- CO2: AEIC 2019 aircraft emissions ---
#AEIC19_DAILY_CO2 kg/m2/s N Y F%y4-%m2-%d2T00:00:00 none none FUELBURN ./HcoDir/AEIC2019/v2022-03/2019/AEIC_2019%m2%d2.0.5x0.625.36L.nc
AEIC19_MONMEAN_CO2 kg/m2/s 2019 Y F%y4-%m2-%d2T00:00:00 none none FUELBURN ./HcoDir/AEIC2019/v2022-03/2019_monmean/AEIC_monmean_2019%m2.0.5x0.625.36L.nc

# --- CO2: Surface correction for CO oxidation ---
FOSSILCO2_MONTHLY kg/m/s N    Y F%y4-%m2-01T01:00:00 none none CO2     ./HcoDir/CO2/v2022-11/FOSSIL/ODIAC_CO2.monthly.generic.1x1.nc
CO2_LIVESTOCK     kg/m/s 2004 Y F%y4-%m2-01T01:00:00 none none CH4_004 ./HcoDir/CO2/v2022-11/CHEM/CH4_source.geos.2x25.nc
CO2_WASTE         kg/m/s 2004 Y F%y4-%m2-01T01:00:00 none none CH4_005 ./HcoDir/CO2/v2022-11/CHEM/CH4_source.geos.2x25.nc
CO2_RICE          kg/m/s 2004 Y F%y4-%m2-01T01:00:00 none none CH4_007 ./HcoDir/CO2/v2022-11/CHEM/CH4_source.geos.2x25.nc
CO2_WETLANDS      kg/m/s 2004 Y F%y4-%m2-01T01:00:00 none none CH4_010 ./HcoDir/CO2/v2022-11/CHEM/CH4_source.geos.2x25.nc
CO2_NATURAL       kg/m/s 2004 Y F%y4-%m2-01T01:00:00 none none CH4_012 ./HcoDir/CO2/v2022-11/CHEM/CH4_source.geos.2x25.nc
CO2_ISOPRENE      kg/m/s 2004 Y F%y4-%m2-01T01:00:00 none none ISOP    ./HcoDir/CO2/v2022-11/CHEM/Isoprene-2004.geos.2x25.nc
CO2_MONOTERP      kg/m/s 2004 Y F%y4-%m2-01T01:00:00 none none MONOT   ./HcoDir/CO2/v2022-11/CHEM/Monoterpene-2004.geos.2x25.nc

# --- CO2 production from CO ---
#CO2_COPROD        kgC/m3/s N Y F%y2-%m2-01T00:00:00 none none LCO  ./HcoDir/CO2/v2022-11/CHEM/CO2_prod_rates.GEOS5.2x25.47L.nc

#==============================================================================
# --- OCS flux data ----
#==============================================================================
OCS_ANTHRO        kgS/km2/s N Y F%y4-%m2-01T00:00:00 none none COS_Flux ./HcoDir/OCS_FLUX/v2022-11/anthro_v3/v3_anthro/%y4/%m2.nc
OCS_BIOMASS       kgS/km2/s N Y F%y4-%m2-01T00:00:00 none none COS_Flux ./HcoDir/OCS_FLUX/v2022-11/Biomassburn/CMSV2-CO2-GFED3-ZC/%y4/%m2.nc
OCS_MISSING_OCEAN kgS/km2/s N Y F%y4-%m2-01T00:00:00 none none COS_Flux ./HcoDir/OCS_FLUX/v2022-11/MissingOcean/LUKAI/%y4/%m2.nc
OCS_OCEAN         kgS/km2/s N Y F%y4-%m2-01T00:00:00 none none COS_Flux ./HcoDir/OCS_FLUX/v2022-11/OceanCOS/Kettle/%y4/%m2.nc
#
###############################################################################
###
### HEMCO Extension Data (update if HEMCO_Config.rc changes)
###
###############################################################################

# --- GFED4 (Extension 111)---
GFED_TEMP kgDM/m2/s N Y F%y4-%m2-01T00:00:00 none none DM_TEMP ./HcoDir/GFED4/v2022-11/%y4/GFED4_gen.025x025.%y4%m2.nc
GFED_AGRI kgDM/m2/s N Y F%y4-%m2-01T00:00:00 none none DM_AGRI ./HcoDir/GFED4/v2022-11/%y4/GFED4_gen.025x025.%y4%m2.nc
GFED_DEFO kgDM/m2/s N Y F%y4-%m2-01T00:00:00 none none DM_DEFO ./HcoDir/GFED4/v2022-11/%y4/GFED4_gen.025x025.%y4%m2.nc
GFED_BORF kgDM/m2/s N Y F%y4-%m2-01T00:00:00 none none DM_BORF ./HcoDir/GFED4/v2022-11/%y4/GFED4_gen.025x025.%y4%m2.nc
GFED_PEAT kgDM/m2/s N Y F%y4-%m2-01T00:00:00 none none DM_PEAT ./HcoDir/GFED4/v2022-11/%y4/GFED4_gen.025x025.%y4%m2.nc
GFED_SAVA kgDM/m2/s N Y F%y4-%m2-01T00:00:00 none none DM_SAVA ./HcoDir/GFED4/v2022-11/%y4/GFED4_gen.025x025.%y4%m2.nc

# --- GFED_daily scale factors ---
GFED_FRAC_DAY 1 N Y %y4-%m2-%d2T00:00:00 none none GFED_FRACDAY ./HcoDir/GFED4/v2022-11/%y4/GFED4_dailyfrac_gen.025x025.%y4%m2.nc

# --- GFED_3hourly (off by default in HEMCO_Config.rc) ---
#GFED_FRAC_3HOUR 1 N Y %y4-%m2-01T%h2:00:00 none none GFED_FRAC3HR ./HcoDir/GFED4/v2020-02/%y4/GFED4_3hrfrac_gen.025x025.%y4%m2.nc

###############################################################################
###
### HEMCO Non-Emissions Data (update if HEMCO_Config.rc changes)
###
###############################################################################

# --- Time zones (offset to UTC) ---
TIMEZONES count N V - none none UTC_OFFSET ./HcoDir/TIMEZONES/v2015-02/timezones_voronoi_1x1.nc

# --- OH from the latest 10-year benchmark ---
# ** TODO: Update to 14.0.0 **
#GLOBAL_OH kg/m3 N N %y4-%m2-01T00:00:00       none none SpeciesConc_OH ./HcoDir/GCClassic_Output/13.0.0/%y4/GEOSChem.SpeciesConc.%y4%m2%d2_0000z.nc4

# --- OH from GEOS-Chem v5-07-08 ---
# --- NOTE: Use this for CH4/IMI ---
#GLOBAL_OH kg/m3 1985 Y F%y4-%m2-01T00:00:00 none none OH ./HcoDir/OH/v2022-11/v5-07-08/OH_3Dglobal.geos5.47L.4x5.nc

# --- Global Cl concentrations
GLOBAL_Cl mol/mol N Y F%y4-%m2-%d2T00:00:00  none none SpeciesConc_Cl ./HcoDir/GCClassic_Output/13.0.0/%y4/GEOSChem.SpeciesConc.%y4%m201_0000z.nc4

# --- Stratospheric L(CO) from GMI ---
CH4_LOSS s-1 1985 Y F%y4-%m2-01T00:00:00 none none CH4loss ./HcoDir/CH4/v2022-11/4x5/gmi.ch4loss.geos5_47L.4x5.nc

# --- P(CO) from CH4 and NMVOC from the last 10-yr benchmark ---
PCO_CH4    molec/cm3/s N Y F%y4-%m2-01T00:00:00 none none ProdCOfromCH4   ./HcoDir/GCClassic_Output/13.0.0/%y4/GEOSChem.ProdLoss.%y4%m201_0000z.nc4
PCO_NMVOC  molec/cm3/s N Y F%y4-%m2-01T00:00:00 none none ProdCOfromNMVOC ./HcoDir/GCClassic_Output/13.0.0/%y4/GEOSChem.ProdLoss.%y4%m201_0000z.nc4

#--- GMI chemistry: prod/loss rates (GMI_PROD_LOSS) ---
GMI_LOSS_CH4 s-1   2000 Y F%y4-%m2-01T00:00:00 none none loss ./HcoDir/GMI/v2022-11/gmi.clim.CH4.geos5.2x25.nc
GMI_PROD_CH4 v/v/s 2000 Y F%y4-%m2-01T00:00:00 none none prod ./HcoDir/GMI/v2022-11/gmi.clim.CH4.geos5.2x25.nc
GMI_LOSS_CO  s-1   2005 Y F%y4-%m2-01T00:00:00 none none loss ./HcoDir/GMI/v2022-11/gmi.clim.CO.geos5.2x25.nc
GMI_PROD_CO  v/v/s 2005 Y F%y4-%m2-01T00:00:00 none none prod ./HcoDir/GMI/v2022-11/gmi.clim.CO.geos5.2x25.nc

###############################################################################
###
### HEMCO Scale Factors (update if HEMCO_Config.rc changes)
###
###############################################################################

# --- CO2 diurnal and weekly scale factors ---
CO2_DIURNAL 1 N    Y F%y4-%m2-%d2T%h2:00:00 none none diurnal_scale_factors ./HcoDir/CO2/v2015-04/FOSSIL/TIMES_diurnal_scale_factors.nc
CO2_WEEKLY  1 2006 Y F%y4-%m2-%d2T00:00:00  none none weekly_scale_factors  ./HcoDir/CO2/v2015-04/FOSSIL/TIMES_weekly_scale_factors.nc

# --- CH4 manure and rice scale factors ---
MANURE_SF  1 N    Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2017-10/Seasonal_SF/EMICH4_Manure_ScalingFactors.WithClimatology.nc
RICE_SF    1 2012 Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2017-10/Seasonal_SF/EMICH4_Rice_ScalingFactors.SetMissing.nc

#--- annual scale factors ---
## Need DC0360xPC0181_CFnnnnx6C.bin
#LIQFUEL_THISYR    1 N Y F%y4-01-01T00:00:00  none none COscalar  ./HcoDir/AnnualScalar/v2014-07/AnnualScalar.geos.1x1.nc
#LIQFUEL_1985      1 N Y F1985-01-01T00:00:00 none none COscalar  ./HcoDir/AnnualScalar/v2014-07/AnnualScalar.geos.1x1.nc
#LIQFUEL_2006      1 N Y F2006-01-01T00:00:00 none none COscalar  ./HcoDir/AnnualScalar/v2014-07/AnnualScalar.geos.1x1.nc
#LIQFUEL_2002      1 N Y F2002-01-01T00:00:00 none none COscalar  ./HcoDir/AnnualScalar/v2014-07/AnnualScalar.geos.1x1.nc
#LIQFUEL_2008_2010 1 N Y F%y4-01-01T00:00:00  none none COscalar  ./HcoDir/AnnualScalar/v2014-07/AnnualScalar.geos.1x1.nc

# --- CEDS vertical partitioning ---
ENERGY_LEVS   1 N Y - none none g_energy   ./HcoDir/VerticalScaleFactors/v2021-05/gc_layers.nc
INDUSTRY_LEVS 1 N Y - none none g_industry ./HcoDir/VerticalScaleFactors/v2021-05/gc_layers.nc
SHIP_LEVS     1 N Y - none none cmv_c3     ./HcoDir/VerticalScaleFactors/v2021-05/gc_layers.nc

# --- CO2 scale factors
AVIATION_SURF_CORR      1       2004 Y F%y2-01-01T00:00:00  none none CO2 ./HcoDir/CO2/v2022-11/FOSSIL/Aviation_SurfCorr_SclFac.1x1.nc
OCEANCO2_SCALED_MONTHLY kg/m2/s N    Y F%y2-%m2-01T00:00:00 none none CO2 ./HcoDir/CO2/v2022-11/OCEAN/Scaled_Ocean_CO2_monthly.nc

###############################################################################
###
### Masks (update if HEMCO_Config.rc changes)
###
###############################################################################

CAN_MASK          1 N Y - none none MASK ./HcoDir/MASKS/v2018-09/Canada_Mask.001x001.nc
CAN_MASK_MIRROR   1 N Y - none none MASK ./HcoDir/MASKS/v2018-09/Canada_Mask_Mirror.001x001.nc
CHINA_MASK        1 N Y - none none MASK ./HcoDir/MASKS/v2018-09/China_mask.generic.1x1.nc
CONUS_MASK        1 N Y - none none MASK ./HcoDir/MASKS/v2018-09/CONUS_Mask.001x001.nc
CONUS_MASK_MIRROR 1 N Y - none none MASK ./HcoDir/MASKS/v2018-09/CONUS_Mask_Mirror.001x001.nc
MEX_MASK          1 N Y - none none MASK ./HcoDir/MASKS/v2018-09/Mexico_Mask.001x001.nc
MEX_MASK_MIRROR   1 N Y - none none MASK ./HcoDir/MASKS/v2018-09/Mexico_Mask_Mirror.001x001.nc
#
%%

DerivedExports%%
# ---------|---------|--------------------------------------------|
#  Export  | Primary |_________________ Mask _____________________|
#  Name    |  Name   |    Name    |           Expression          |
# ---------|---------|------------|-------------------------------|
# ---------|---------|------------|-------------------------------|
${RUNDIR_MET_EXTDATA_DERIVED_EXPORTS}
%%
