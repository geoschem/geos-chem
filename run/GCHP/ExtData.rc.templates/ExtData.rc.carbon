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
#
# Comment all entries in this section when
# "--> USE_CH4_DATA : false" is in HEMCO_Config.rc.
#==============================================================================

# --- CH4: Gridded GHGI v2 Express Extension (Maasakkers et al., submitted to ES&T, 2023) ---
#
# Comment these lines when "GHGI_v2_Express_Ext: false" is in HEMCO_Config.rc.
GHGI_EE_OIL_EXPLORATION         molec/cm2/s N    Y - none none  emi_ch4_1B2a_Petroleum_Systems_Exploration   ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_OIL_PRODUCTION          molec/cm2/s N    Y - none none  emi_ch4_1B2a_Petroleum_Systems_Production    ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_OIL_REFINING            molec/cm2/s N    Y - none none  emi_ch4_1B2a_Petroleum_Systems_Refining      ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_OIL_TRANSPORT           molec/cm2/s N    Y - none none  emi_ch4_1B2a_Petroleum_Systems_Transport     ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_GAS_DISTRIBUTION        molec/cm2/s N    Y - none none  emi_ch4_1B2b_Natural_Gas_Distribution        ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_GAS_EXPLORATION         molec/cm2/s N    Y - none none  emi_ch4_1B2b_Natural_Gas_Exploration         ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_GAS_PROCESSING          molec/cm2/s N    Y - none none  emi_ch4_1B2b_Natural_Gas_Processing          ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_GAS_PRODUCTION          molec/cm2/s N    Y - none none  emi_ch4_1B2b_Natural_Gas_Production          ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_GAS_TRANSMISSION        molec/cm2/s N    Y - none none  emi_ch4_1B2b_Natural_Gas_TransmissionStorage ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAL_UNDERGROUND        molec/cm2/s N    Y - none none  emi_ch4_1B1a_Underground_Coal                ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAL_SURFACE            molec/cm2/s N    Y - none none  emi_ch4_1B1a_Surface_Coal                    ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAL_ABANDONED          molec/cm2/s N    Y - none none  emi_ch4_1B1a_Abandoned_Coal                  ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_LIVESTOCK_ENT           molec/cm2/s N    Y - none none  emi_ch4_3A_Enteric_Fermentation              ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_LIVESTOCK_MAN           molec/cm2/s N    Y - none none  emi_ch4_3B_Manure_Management                 ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_LANDFILLS_IND           molec/cm2/s N    Y - none none  emi_ch4_5A1_Landfills_Industrial             ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_LANDFILLS_MSW           molec/cm2/s N    Y - none none  emi_ch4_5A1_Landfills_MSW                    ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_LANDFILLS_COMP          molec/cm2/s N    Y - none none  emi_ch4_5B1_Composting                       ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_WASTEWATER_DOM          molec/cm2/s N    Y - none none  emi_ch4_5D_Wastewater_Treatment_Domestic     ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_WASTEWATER_IND          molec/cm2/s N    Y - none none  emi_ch4_5D_Wastewater_Treatment_Industrial   ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_RICE                    molec/cm2/s N    Y - none none  emi_ch4_3C_Rice_Cultivation                  ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_OTHER__MCOMB            molec/cm2/s N    Y - none none  emi_ch4_1A_Combustion_Mobile                 ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_OTHER__SCOMB            molec/cm2/s N    Y - none none  emi_ch4_1A_Combustion_Stationary             ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_OTHER__PIND             molec/cm2/s N    Y - none none  emi_ch4_2B8_Industry_Petrochemical           ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_OTHER__FIND             molec/cm2/s N    Y - none none  emi_ch4_2C2_Industry_Ferroalloy              ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_OTHER__BURN             molec/cm2/s N    Y - none none  emi_ch4_3F_Field_Burning                     ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_OTHER__ABOG             molec/cm2/s N    Y - none none  emi_ch4_1B2ab_Abandoned_Oil_Gas              ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAST_OIL_EXPLORATION   molec/cm2/s N    Y - none none  emi_ch4_1B2a_Petroleum_Systems_Exploration   ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAST_OIL_PRODUCTION    molec/cm2/s N    Y - none none  emi_ch4_1B2a_Petroleum_Systems_Production    ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAST_OIL_REFINING      molec/cm2/s N    Y - none none  emi_ch4_1B2a_Petroleum_Systems_Refining      ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAST_OIL_TRANSPORT     molec/cm2/s N    Y - none none  emi_ch4_1B2a_Petroleum_Systems_Transport     ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAST_GAS_DISTRIBUTION  molec/cm2/s N    Y - none none  emi_ch4_1B2b_Natural_Gas_Distribution        ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAST_GAS_EXPLORATION   molec/cm2/s N    Y - none none  emi_ch4_1B2b_Natural_Gas_Exploration         ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAST_GAS_PROCESSING    molec/cm2/s N    Y - none none  emi_ch4_1B2b_Natural_Gas_Processing          ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAST_GAS_PRODUCTION    molec/cm2/s N    Y - none none  emi_ch4_1B2b_Natural_Gas_Production          ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAST_GAS_TRANSMISSION  molec/cm2/s N    Y - none none  emi_ch4_1B2b_Natural_Gas_TransmissionStorage ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAST_COAL_UNDERGROUND  molec/cm2/s N    Y - none none  emi_ch4_1B1a_Underground_Coal                ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAST_COAL_SURFACE      molec/cm2/s N    Y - none none  emi_ch4_1B1a_Surface_Coal                    ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAST_COAL_ABANDONED    molec/cm2/s N    Y - none none  emi_ch4_1B1a_Abandoned_Coal                  ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAST_LIVESTOCK_ENT     molec/cm2/s N    Y - none none  emi_ch4_3A_Enteric_Fermentation              ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAST_LIVESTOCK_MAN     molec/cm2/s N    Y - none none  emi_ch4_3B_Manure_Management                 ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAST_LANDFILLS_IND     molec/cm2/s N    Y - none none  emi_ch4_5A1_Landfills_Industrial             ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAST_LANDFILLS_MSW     molec/cm2/s N    Y - none none  emi_ch4_5A1_Landfills_MSW                    ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAST_LANDFILLS_COMP    molec/cm2/s N    Y - none none  emi_ch4_5B1_Composting                       ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAST_WASTEWATER_DOM    molec/cm2/s N    Y - none none  emi_ch4_5D_Wastewater_Treatment_Domestic     ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAST_WASTEWATER_IND    molec/cm2/s N    Y - none none  emi_ch4_5D_Wastewater_Treatment_Industrial   ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAST_RICE              molec/cm2/s N    Y - none none  emi_ch4_3C_Rice_Cultivation                  ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAST_OTHER__MCOMB      molec/cm2/s N    Y - none none  emi_ch4_1A_Combustion_Mobile                 ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAST_OTHER__SCOMB      molec/cm2/s N    Y - none none  emi_ch4_1A_Combustion_Stationary             ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAST_OTHER__PIND       molec/cm2/s N    Y - none none  emi_ch4_2B8_Industry_Petrochemical           ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAST_OTHER__FIND       molec/cm2/s N    Y - none none  emi_ch4_2C2_Industry_Ferroalloy              ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAST_OTHER__BURN       molec/cm2/s N    Y - none none  emi_ch4_3F_Field_Burning                     ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAST_OTHER__ABOG       molec/cm2/s N    Y - none none  emi_ch4_1B2ab_Abandoned_Oil_Gas              ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc


# --- Scarpelli Mexico ---
#
# Comment these when "--> Scarpelli_Mexico : false" is in HEMCO_Config.rc.
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
#
# Comment these when "--> Scarpelli_Canada : false" is in HEMCO_Config.rc
#
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
#
# Comment these when "--> GFEIv2 : false" is in HEMCO_Config.rc
GFEI_CH4_OIL  kg/m2/s N Y - none none  emis_ch4  ./HcoDir/CH4/v2022-01/GFEIv2/Global_Fuel_Exploitation_Inventory_v2_2019_Oil_All.nc
GFEI_CH4_GAS  kg/m2/s N Y - none none  emis_ch4  ./HcoDir/CH4/v2022-01/GFEIv2/Global_Fuel_Exploitation_Inventory_v2_2019_Gas_All.nc
GFEI_CH4_COAL kg/m2/s N Y - none none  emis_ch4  ./HcoDir/CH4/v2022-01/GFEIv2/Global_Fuel_Exploitation_Inventory_v2_2019_Coal.nc


# --- EDGAR v7.0 ---
#
# Use updated files (v2024-01); these are COARDS-compliant.
#
# Comment these when "--> EDGARv7 : false" is in HEMCO_Config.rc
EDGAR7_CH4_OILGASCOAL__1B1_1B2  kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2024-01/EDGARv7/%y4/v7.0_FT2021_CH4_%y4_PRO.0.1x0.1.nc
EDGAR7_CH4_OTHER__1A1_1B1_1B2   kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2024-01/EDGARv7/%y4/v7.0_FT2021_CH4_%y4_REF_TRF.0.1x0.1.nc
EDGAR7_CH4_LIVESTOCK__4A        kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2024-01/EDGARv7/%y4/v7.0_FT2021_CH4_%y4_ENF.0.1x0.1.nc
EDGAR7_CH4_LIVESTOCK__4B        kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2024-01/EDGARv7/%y4/v7.0_FT2021_CH4_%y4_MNM.0.1x0.1.nc
EDGAR7_CH4_LANDFILLS__6A_6D     kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2024-01/EDGARv7/%y4/v7.0_FT2021_CH4_%y4_SWD_LDF.0.1x0.1.nc
EDGAR7_CH4_WASTEWATER__6B       kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2024-01/EDGARv7/%y4/v7.0_FT2021_CH4_%y4_WWT.0.1x0.1.nc
EDGAR7_CH4_RICE__4C_4D          kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2024-01/EDGARv7/%y4/v7.0_FT2021_CH4_%y4_AGS.0.1x0.1.nc
EDGAR7_CH4_OTHER__1A1a          kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2024-01/EDGARv7/%y4/v7.0_FT2021_CH4_%y4_ENE.0.1x0.1.nc
EDGAR7_CH4_OTHER__1A2           kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2024-01/EDGARv7/%y4/v7.0_FT2021_CH4_%y4_IND.0.1x0.1.nc
EDGAR7_CH4_OTHER__1A3a_CDS      kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2024-01/EDGARv7/%y4/v7.0_FT2021_CH4_%y4_TNR_Aviation_CDS.0.1x0.1.nc
EDGAR7_CH4_OTHER__1A3a_CRS      kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2024-01/EDGARv7/%y4/v7.0_FT2021_CH4_%y4_TNR_Aviation_CRS.0.1x0.1.nc
EDGAR7_CH4_OTHER__1A3a_LTO      kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2024-01/EDGARv7/%y4/v7.0_FT2021_CH4_%y4_TNR_Aviation_LTO.0.1x0.1.nc
EDGAR7_CH4_OTHER__1A3b          kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2024-01/EDGARv7/%y4/v7.0_FT2021_CH4_%y4_TRO_noRES.0.1x0.1.nc
EDGAR7_CH4_OTHER__1A3c_1A3e     kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2024-01/EDGARv7/%y4/v7.0_FT2021_CH4_%y4_TNR_Other.0.1x0.1.nc
EDGAR7_CH4_OTHER__1A3d_1C2      kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2024-01/EDGARv7/%y4/v7.0_FT2021_CH4_%y4_TNR_Ship.0.1x0.1.nc
EDGAR7_CH4_OTHER__1A4           kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2024-01/EDGARv7/%y4/v7.0_FT2021_CH4_%y4_RCO.0.1x0.1.nc
EDGAR7_CH4_OTHER__2B            kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2024-01/EDGARv7/%y4/v7.0_FT2021_CH4_%y4_CHE.0.1x0.1.nc
EDGAR7_CH4_OTHER__2C            kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2024-01/EDGARv7/%y4/v7.0_FT2021_CH4_%y4_IRO.0.1x0.1.nc
EDGAR7_CH4_OTHER__4F            kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2024-01/EDGARv7/%y4/v7.0_FT2021_CH4_%y4_AWB.0.1x0.1.nc
EDGAR7_CH4_OTHER__6C            kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2024-01/EDGARv7/%y4/v7.0_FT2021_CH4_%y4_SWD_INC.0.1x0.1.nc


# --- JPL WetCHARTs v1.0 (Bloom et al., https://doi.org/10.3334/ORNLDAAC/1502) ---
#
# Use updated files (v2024-01); these are COARDS-compliant.
#
# Comment this when "--> JPL_WETCHARTS: false" is in HEMCO_Config.rc.
JPLW_CH4  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none  emi_ch4  ./HcoDir/CH4/v2024-01/JPL_WetCharts/HEensemble/JPL_WetCharts_2010-2019.Ensemble_Mean.0.5x0.5.nc


# --- Geological seeps ---
#
# Comment this when "--> SEEPS : false" is in HEMCO_Config.rc.
CH4_SEEPS kg/m2/s N Y - none none  emi_ch4  ./HcoDir/CH4/v2020-04/Seeps/Etiope_CH4GeologicalEmis_ScaledToHmiel.1x1.nc


# --- Emissions from reservoirs (Delwich et al 2022) ---
#
# Use updated files (v2024-01); these are COARDS-compliant.
#
# Comment these when "--> RESERVOIRS : false" is in HEMCO_Config.rc.
CH4_RES_DAM kg/m2/s N Y - none none CH4emis  ./HcoDir/CH4/v2024-01/ResME/ResME_Dam_Emissions.0.1x0.1.nc
CH4_RES_SFC kg/m2/s N Y - none none CH4emis  ./HcoDir/CH4/v2024-01/ResME/ResME_Surface_Emissions.0.1x0.1.nc


# --- Emissions from termites (Fung et al 1991) ---
#
# Comment this when "--> FUNG_TERMITES : false" is in HEMCO_Config.rc.
CH4_TERMITES    kg/m2/s N Y - none none  CH4  ./HcoDir/CH4/v2022-11/4x5/termites.geos.4x5.nc


# --- Soil absorption fro MeMo model (Murguia-Flores et al. 2018, GMD) ---
#
# Note; Use climatology by default
# Comment this when "--> MeMo_SOIL_ABSORPTION :  false" is in HEMCO_Config.rc.
CH4_SOILABSORB kg/m2/s  2009 Y F%y4-%m2-01T00:00:00 none none  CH4uptake  ./HcoDir/CH4/v2019-10/MeMo_SoilAbs/MeMo_CH4uptake_Climatology.nc
#
# Annual emissions are disabled by default
#CH4_SOILABSORB kg/m2/s N    Y F%y4-%m2-01T00:00:00 none none  CH4uptake  ./HcoDir/CH4/v2019-10/MeMo_SoilAbs/MeMo_CH4uptake_1990-2009.nc

# --- BBMIPs ---
# Add entries here

# --- CEDS (historical) or Shared Socioeconomic Pathways (future) ---
# Add entries here

# --- EDGAR v6.0 ---
# Add entries here

# --- Gridded GHGI emissions ---
# Add entries here

# --- QFED v2.4r8 ---
# Add entries here

# --- Emissions from lakes (Maasakkers et al 2019) ---
# Add entries here

# --- Soil absoprtion (Fung et al 1991) ---
# Add entries here

# --- RCP future emission scenarios ---
# Add entries here


#==============================================================================
# %%%%% CO emissions %%%%%
#
# Comment all entries in this section when
# "--> USE_CO_DATA : false" is in HEMCO_Config.rc.
#==============================================================================

# --- AEIC 2019 aircraft (AEIC) ---
#
# Comment this when "--> AEIC2019_MONMEAN :  false" is in HEMCO_Config.rc.
AEIC19_MONMEAN_CO kg/m2/s Y Y F2019-%m2-01T00:00:00 none none CO ./HcoDir/AEIC2019/v2022-03/2019_monmean/AEIC_monmean_2019%m2.0.5x0.625.36L.nc
#
# Comment this when "--> AEIC2019_DAILY :  false" is in HEMCO_Config.rc.
##AEIC19_DAILY_CO   kg/m2/s 2019 Y F%y4-%m2-%d2T00:00:00 none none CO ./HcoDir/AEIC2019/v2022-03/2019/%m2/AEIC_2019%m2%d2.0.5x0.625.36L.nc


# --- CEDSv2 emissions ---
#
# Comment these when "--> CEDSv2 : false" is in HEMCO_Config.rc.
CEDS_CO_AGR  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO_agr  ./HcoDir/CEDS/v2021-06/%y4/CO-em-anthro_CMIP_CEDS_%y4.nc
CEDS_CO_ENE  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO_ene  ./HcoDir/CEDS/v2021-06/%y4/CO-em-anthro_CMIP_CEDS_%y4.nc
CEDS_CO_IND  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO_ind  ./HcoDir/CEDS/v2021-06/%y4/CO-em-anthro_CMIP_CEDS_%y4.nc
CEDS_CO_TRA  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO_tra  ./HcoDir/CEDS/v2021-06/%y4/CO-em-anthro_CMIP_CEDS_%y4.nc
CEDS_CO_RCO  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO_rco  ./HcoDir/CEDS/v2021-06/%y4/CO-em-anthro_CMIP_CEDS_%y4.nc
CEDS_CO_SLV  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO_slv  ./HcoDir/CEDS/v2021-06/%y4/CO-em-anthro_CMIP_CEDS_%y4.nc
CEDS_CO_WST  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO_wst  ./HcoDir/CEDS/v2021-06/%y4/CO-em-anthro_CMIP_CEDS_%y4.nc
#
# Comment this when  "--> CEDSv2_SHIP : false" is in HEMCO_Config.rc.
CEDS_CO_SHP  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO_shp  ./HcoDir/CEDS/v2021-06/%y4/CO-em-anthro_CMIP_CEDS_%y4.nc


# --- APEI (Canada) ---
# Add entries here

# --- BB4MIPs ---
# Add entries here

# --- CEDS GBD-MAPS ---
# Add entries here

# --- CMIP6 ---
# Add entries here

# --- DICE-Africa ---
# Add entries here

# --- GFAS biomass burning ---
# Add entries here

# --- HTAP ---
# Add entries here

# --- NEI2016_MONMEAN (USA) ---
# Add entries here

# --- QFED2 biomass burning ---
# Add entries here

# --- RCP future emissions scenarios ---
# Add entries here


#==============================================================================
# --- CO2 emissions ---
#
# Comment all entries in this section when
# "--> USE_CO2_DATA : false" is in HEMCO_Config.rc.
#==============================================================================

# --- CO2 Fossil Fuel ---
#
# Comment this when "--> FOSSIL_CDIAC : false" is in HEMCO_Config.rc.
##FOSSILCO2_CDIAC  kg/m2/s N Y  F%y4-%m2-01T00:00:00 none none CO2 ./HcoDir/CO2/v2014-09/FOSSIL/CDIAC_v2014.monthly.generic.1x1.nc
#
# Comment this when "--> FOSSIL_ODIAC : false" is in HEMCO_Config.rc.
FOSSILCO2_ODIAC  kg/m2/s N Y  F%y4-%m2-01T00:00:00 none none CO2 ./HcoDir/CO2/v2022-11/FOSSIL/ODIAC_CO2.monthly.generic.1x1.nc


# --- CO2 Ocean Exchange ---
#
# Comment this when "--> OCEAN_EXCH_TAKA09 : false" is in HEMCO_Config.rc.
##OCEANCO2_TAKA_ANNUAL   kg/m2/s N Y - none none CO2 ./HcoDir/CO2/v2022-11/OCEAN/Taka2009_CO2_Annual.nc
#
# Comment this when "--> OCEAN_EXCH_SCALED : false" is in HEMCO_Config.rc.
OCEANCO2_TAKA_MONTHLY  kg/m2/s 2000 Y F%y4-%m2-01T00:00:00 none none CO2 ./HcoDir/CO2/v2022-11/OCEAN/Taka2009_CO2_Monthly.nc


# --- Balanced biosphere exchange ---
#
# Comment this when "--> BBIO_DIURNAL : false" is in HEMCO_Config.rc.
##BBIOCO2_DIURNAL kg/m2/s N Y F%y4-%m2-%d2T%h2-01-01 none none CO2 ./HcoDir/CO2/v2014-09/BBIO/BBIO_diurnal_CO2.nc
#
# Comment this when "--> BBIO_SIB3 : false" is in HEMCO_Config.rc.
SIB_BBIO_CO2 kg/m2/s N Y F%y4-%m2-%d2T%h2-01-01 none none CO2 ./HcoDir/CO2/v2022-11/BIO/SiB3_3hr_NEP.nc


# --- Net Terrestrial Exchange ---
#
# Comment this when "--> NET_TERR_EXCH : false" is in HEMCO_Config.rc.
CO2_NET_TERRESTRIAL kg/m2/s N Y - none none CO2 ./HcoDir/CO2/v2022-11/BIO/Net_terrestrial_exch_5.29Pg.generic.1x1.nc


# --- CO2 Ship Emissions ---
#
# Comment this when "--> CEDS_v2_SHIP : false" is in HEMCO_Config.rc.
CEDS_CO2_SHP kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO2_shp ./HcoDir/CEDS/v2021-06/%y4/CO2-em-anthro_CMIP_CEDS_%y4.nc


# --- CO2: AEIC 2019 aircraft emissions ---
## Comment this when "--> AEIC2019_DAILY :  false" is in HEMCO_Config.rc.
##AEIC19_DAILY_CO2   kg/m2/s 2019 Y F%y4-%m2-%d2T00:00:00 none none FUELBURN ./HcoDir/AEIC2019/v2022-03/2019/%m2/AEIC_2019%m2%d2.0.5x0.625.36L.nc
#
# Comment this when "--> AEIC2019_MONMEAN :  false" is in HEMCO_Config.rc.
AEIC19_MONMEAN_CO2 kg/m2/s 2019 Y F%y4-%m2-%d2T00:00:00 none none FUELBURN ./HcoDir/AEIC2019/v2022-03/2019_monmean/AEIC_monmean_2019%m2.0.5x0.625.36L.nc


# --- CO2: Surface correction for CO oxidation ---
#
# Comment these when "--> CO2CORR : false" is in HEMCO_Config.rc.
FOSSILCO2_MONTHLY kg/m/s N    Y F%y4-%m2-01T01:00:00 none none CO2     ./HcoDir/CO2/v2022-11/FOSSIL/ODIAC_CO2.monthly.generic.1x1.nc
CO2_LIVESTOCK     kg/m/s 2004 Y F%y4-%m2-01T01:00:00 none none CH4_004 ./HcoDir/CO2/v2022-11/CHEM/CH4_source.geos.2x25.nc
CO2_WASTE         kg/m/s 2004 Y F%y4-%m2-01T01:00:00 none none CH4_005 ./HcoDir/CO2/v2022-11/CHEM/CH4_source.geos.2x25.nc
CO2_RICE          kg/m/s 2004 Y F%y4-%m2-01T01:00:00 none none CH4_007 ./HcoDir/CO2/v2022-11/CHEM/CH4_source.geos.2x25.nc
CO2_WETLANDS      kg/m/s 2004 Y F%y4-%m2-01T01:00:00 none none CH4_010 ./HcoDir/CO2/v2022-11/CHEM/CH4_source.geos.2x25.nc
CO2_NATURAL       kg/m/s 2004 Y F%y4-%m2-01T01:00:00 none none CH4_012 ./HcoDir/CO2/v2022-11/CHEM/CH4_source.geos.2x25.nc
CO2_ISOPRENE      kg/m/s 2004 Y F%y4-%m2-01T01:00:00 none none ISOP    ./HcoDir/CO2/v2022-11/CHEM/Isoprene-2004.geos.2x25.nc
CO2_MONOTERP      kg/m/s 2004 Y F%y4-%m2-01T01:00:00 none none MONOT   ./HcoDir/CO2/v2022-11/CHEM/Monoterpene-2004.geos.2x25.nc


# --- CO2 production from CO ---
#
# Use updated files (v2024-01); these are COARDS-compliant.
#
# Comment this when "--> CO2_COPROD : false" is in HEMCO_Config.rc.
CO2_COPROD        kgC/m3/s N Y F%y2-%m2-01T00:00:00 none none LCO  ./HcoDir/CO2/v2024-01/CHEM/CO2_prod_rates.GEOS5.2x25.72L.nc

#==============================================================================
# --- OCS emission fluxes ---
#
# If USE_OCS_DATA is set to "false" in HEMCO_Config.rc,
# then you must comment out all of these entries.
#==============================================================================
OCS_ANTHRO        kgS/km2/s N Y F%y4-%m2-01T00:00:00 none none COS_Flux ./HcoDir/OCS_FLUX/v2022-11/anthro_v3/v3_anthro/%y4/%m2.nc
OCS_BIOMASS       kgS/km2/s N Y F%y4-%m2-01T00:00:00 none none COS_Flux ./HcoDir/OCS_FLUX/v2022-11/Biomassburn/CMSV2-CO2-GFED3-ZC/%y4/%m2.nc
OCS_MISSING_OCEAN kgS/km2/s N Y F%y4-%m2-01T00:00:00 none none COS_Flux ./HcoDir/OCS_FLUX/v2022-11/MissingOcean/LUKAI/%y4/%m2.nc
OCS_OCEAN         kgS/km2/s N Y F%y4-%m2-01T00:00:00 none none COS_Flux ./HcoDir/OCS_FLUX/v2022-11/OceanCOS/Kettle/%y4/%m2.nc

###############################################################################
###
### HEMCO Extension Data (update if HEMCO_Config.rc changes)
###
###############################################################################

# --- GFED4 (Extension 111)---
#
# Comment these when  "111 GFEE : off" in HEMCO_Config.rc.
GFED_TEMP kgDM/m2/s N Y F%y4-%m2-01T00:00:00 none none DM_TEMP ./HcoDir/GFED4/v2023-03/%y4/GFED4_gen.025x025.%y4%m2.nc
GFED_AGRI kgDM/m2/s N Y F%y4-%m2-01T00:00:00 none none DM_AGRI ./HcoDir/GFED4/v2023-03/%y4/GFED4_gen.025x025.%y4%m2.nc
GFED_DEFO kgDM/m2/s N Y F%y4-%m2-01T00:00:00 none none DM_DEFO ./HcoDir/GFED4/v2023-03/%y4/GFED4_gen.025x025.%y4%m2.nc
GFED_BORF kgDM/m2/s N Y F%y4-%m2-01T00:00:00 none none DM_BORF ./HcoDir/GFED4/v2023-03/%y4/GFED4_gen.025x025.%y4%m2.nc
GFED_PEAT kgDM/m2/s N Y F%y4-%m2-01T00:00:00 none none DM_PEAT ./HcoDir/GFED4/v2023-03/%y4/GFED4_gen.025x025.%y4%m2.nc
GFED_SAVA kgDM/m2/s N Y F%y4-%m2-01T00:00:00 none none DM_SAVA ./HcoDir/GFED4/v2023-03/%y4/GFED4_gen.025x025.%y4%m2.nc
#
# --- GFED_daily scale factors ---
#
# Comment this when "--> GFED_daily : false" is in HEMCO_Config.rc.
GFED_FRAC_DAY 1 N Y %y4-%m2-%d2T00:00:00 none none GFED_FRACDAY ./HcoDir/GFED4/v2023-03/%y4/GFED4_dailyfrac_gen.025x025.%y4%m2.nc
#
# --- GFED_3hourly --
#
# Comment this when "--> GFED_hourly : false" is in HEMCO_Config.rc.
#GFED_FRAC_3HOUR 1 N Y %y4-%m2-01T%h2:00:00 none none GFED_FRAC3HR ./HcoDir/GFED4/v2023-03/%y4/GFED4_3hrfrac_gen.025x025.%y4%m2.nc


###############################################################################
###
### HEMCO Non-Emissions Data (update if HEMCO_Config.rc changes)
###
###############################################################################

# --- Time zones (offset to UTC) ---
TIMEZONES count N V - none none UTC_OFFSET ./HcoDir/TIMEZONES/v2015-02/timezones_voronoi_1x1.nc


# --- OH from the latest 10-year benchmark ---
#
# Comment this if "--> GLOBAL_OH_GC14 : false" is in HEMCO_Config.rc
##GLOBAL_OH kg/m3 N N %y4-%m2-01T00:00:00       none none SpeciesConc_OH ./HcoDir/GCClassic_Output/14.0.0/%y4/GEOSChem.SpeciesConc.%y4%m201_0000z.nc4


# --- OH from GEOS-Chem v5-07-08 ---
# --- NOTE: Use this for CH4/IMI ---
#
# Comment this if "--> GLOBAL_OH_GCv5 : false" is in HEMCO_Config.rc
GLOBAL_OH kg/m3 1985 Y F%y4-%m2-01T00:00:00 none none OH ./HcoDir/OH/v2022-11/v5-07-08/OH_3Dglobal.geos5.72L.4x5.nc


# --- Global Cl concentrations ---
#
# Use updated files (v2024-01); these are COARDS-compliant.
GLOBAL_Cl mol/mol N Y F%y4-%m2-01T00:00:00  none none SpeciesConc_Cl ./HcoDir/CH4/v2024-01/GCC_14_Output/%y4/GEOSChem.SpeciesConc.%y4%m201_0000z.nc4


# --- Stratospheric L(CO) from GMI ---
#
# Use updated files (v2024-01); these are COARDS-compliant.
CH4_LOSS s-1 1985 Y F%y4-%m2-01T00:00:00 none none CH4loss ./HcoDir/CH4/v2024-01/GC_CH4_LOSS/GCC14_72LM.ch4loss.4x5.nc4

# --- P(CO) from CH4 and NMVOC from the last 10-yr benchmark ---
#
# Use updated files (v2024-01); these are COARDS-compliant.
PCO_CH4    molec/cm3/s N Y F%y4-%m2-01T00:00:00 none none ProdCOfromCH4   ./HcoDir/CH4/v2024-01/GCC_14_Output/%y4/GEOSChem.ProdLoss.%y4%m201_0000z.nc4
PCO_NMVOC  molec/cm3/s N Y F%y4-%m2-01T00:00:00 none none ProdCOfromNMVOC ./HcoDir/CH4/v2024-01/GCC_14_Output/%y4/GEOSChem.ProdLoss.%y4%m201_0000z.nc4

#--- GMI chemistry: prod/loss rates (GMI_PROD_LOSS) ---
GMI_LOSS_CO  s-1   2005 Y F%y4-%m2-01T00:00:00 none none loss ./HcoDir/GMI/v2022-11/gmi.clim.CO.geos5.2x25.nc
GMI_PROD_CO  v/v/s 2005 Y F%y4-%m2-01T00:00:00 none none prod ./HcoDir/GMI/v2022-11/gmi.clim.CO.geos5.2x25.nc

###############################################################################
###
### HEMCO Scale Factors (update if HEMCO_Config.rc changes)
###
###############################################################################

#==============================================================================
# %%%%% Scale factors for CH4 emissions %%%%%
#
# Comment out all entries in this section if the USE_CH4_DATA switch
# is set to "false" in HEMCO_Config.rc.
#==============================================================================

# --- GHGI seasonal scale factors ---
#
# Comment these if "--> GHGI_v2_Express_Ext : false" is in HEMCO_Config.rc,
GHGI_OTH_STA_SF  1 N Y F%y4-%m2-%d2T00:00:00 none none monthly_scale_factor_1A_Combustion_Stationary           ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Gridded_GHGI_Methane_v2_Monthly_Scale_Factors_%y4.nc
GHGI_OIL_EXP_SF  1 N Y F%y4-%m2-%d2T00:00:00 none none monthly_scale_factor_1B2a_Petroleum_Systems_Exploration ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Gridded_GHGI_Methane_v2_Monthly_Scale_Factors_%y4.nc
GHGI_OIL_PRD_SF  1 N Y F%y4-%m2-%d2T00:00:00 none none monthly_scale_factor_1B2a_Petroleum_Systems_Production  ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Gridded_GHGI_Methane_v2_Monthly_Scale_Factors_%y4.nc
GHGI_OIL_REF_SF  1 N Y F%y4-%m2-%d2T00:00:00 none none monthly_scale_factor_1B2a_Petroleum_Systems_Refining    ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Gridded_GHGI_Methane_v2_Monthly_Scale_Factors_%y4.nc
GHGI_OIL_TRA_SF  1 N Y F%y4-%m2-%d2T00:00:00 none none monthly_scale_factor_1B2a_Petroleum_Systems_Transport   ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Gridded_GHGI_Methane_v2_Monthly_Scale_Factors_%y4.nc 
GHGI_GAS_EXP_SF  1 N Y F%y4-%m2-%d2T00:00:00 none none monthly_scale_factor_1B2b_Natural_Gas_Exploration       ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Gridded_GHGI_Methane_v2_Monthly_Scale_Factors_%y4.nc 
GHGI_GAS_PRD_SF  1 N Y F%y4-%m2-%d2T00:00:00 none none monthly_scale_factor_1B2b_Natural_Gas_Production        ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Gridded_GHGI_Methane_v2_Monthly_Scale_Factors_%y4.nc
GHGI_LIV_MAN_SF  1 N Y F%y4-%m2-%d2T00:00:00 none none monthly_scale_factor_3B_Manure_Management               ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Gridded_GHGI_Methane_v2_Monthly_Scale_Factors_%y4.nc 
GHGI_RIC_CUL_SF  1 N Y F%y4-%m2-%d2T00:00:00 none none monthly_scale_factor_3C_Rice_Cultivation                ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Gridded_GHGI_Methane_v2_Monthly_Scale_Factors_%y4.nc
GHGI_OTH_BUR_SF  1 N Y F%y4-%m2-%d2T00:00:00 none none monthly_scale_factor_3F_Field_Burning                   ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Gridded_GHGI_Methane_v2_Monthly_Scale_Factors_%y4.nc


# --- Scarpelli_Mexico manure and rice scale factors ---
#
# Comment these if "--> Scarpelli_Mexico : false" is in HEMCO_Config.rc,
MANURE_SF  1 N    Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2017-10/Seasonal_SF/EMICH4_Manure_ScalingFactors.WithClimatology.nc
RICE_SF    1 2012 Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2017-10/Seasonal_SF/EMICH4_Rice_ScalingFactors.SetMissing.nc


# --- Monthly scale factors for EDGARv7 ---
#
# Comment these if "--> EDGARv7 : false" is in HEMCO_Config.rc,
EDGAR_SEASONAL_SF_AGS        1 2018 Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2023-04/EDGARv6_SF/EDGARv6_CH4_MonthlyScaleFactors_AGS.0.1x0.1.nc
EDGAR_SEASONAL_SF_AWB        1 2018 Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2023-04/EDGARv6_SF/EDGARv6_CH4_MonthlyScaleFactors_AWB.0.1x0.1.nc
EDGAR_SEASONAL_SF_CHE        1 2018 Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2023-04/EDGARv6_SF/EDGARv6_CH4_MonthlyScaleFactors_CHE.0.1x0.1.nc
EDGAR_SEASONAL_SF_ENE        1 2018 Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2023-04/EDGARv6_SF/EDGARv6_CH4_MonthlyScaleFactors_ENE.0.1x0.1.nc
EDGAR_SEASONAL_SF_ENF        1 2018 Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2023-04/EDGARv6_SF/EDGARv6_CH4_MonthlyScaleFactors_ENF.0.1x0.1.nc
EDGAR_SEASONAL_SF_FFF        1 2018 Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2023-04/EDGARv6_SF/EDGARv6_CH4_MonthlyScaleFactors_FFF.0.1x0.1.nc
EDGAR_SEASONAL_SF_IND        1 2018 Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2023-04/EDGARv6_SF/EDGARv6_CH4_MonthlyScaleFactors_IND.0.1x0.1.nc
EDGAR_SEASONAL_SF_IRO        1 2018 Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2023-04/EDGARv6_SF/EDGARv6_CH4_MonthlyScaleFactors_IRO.0.1x0.1.nc
EDGAR_SEASONAL_SF_MNM        1 2018 Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2023-04/EDGARv6_SF/EDGARv6_CH4_MonthlyScaleFactors_MNM.0.1x0.1.nc
EDGAR_SEASONAL_SF_PRO        1 2018 Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2023-04/EDGARv6_SF/EDGARv6_CH4_MonthlyScaleFactors_PRO.0.1x0.1.nc
EDGAR_SEASONAL_SF_RCO        1 2018 Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2023-04/EDGARv6_SF/EDGARv6_CH4_MonthlyScaleFactors_RCO.0.1x0.1.nc
EDGAR_SEASONAL_SF_REF_TRF    1 2018 Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2023-04/EDGARv6_SF/EDGARv6_CH4_MonthlyScaleFactors_REF_TRF.0.1x0.1.nc
EDGAR_SEASONAL_SF_SWD_INC    1 2018 Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2023-04/EDGARv6_SF/EDGARv6_CH4_MonthlyScaleFactors_SWD_INC.0.1x0.1.nc
EDGAR_SEASONAL_SF_SWD_LDF    1 2018 Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2023-04/EDGARv6_SF/EDGARv6_CH4_MonthlyScaleFactors_SWD_LDF.0.1x0.1.nc
EDGAR_SEASONAL_SF_TNR_AV_CDS 1 2018 Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2023-04/EDGARv6_SF/EDGARv6_CH4_MonthlyScaleFactors_TNR_Aviation_CDS.0.1x0.1.nc
EDGAR_SEASONAL_SF_TNR_AV_CRS 1 2018 Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2023-04/EDGARv6_SF/EDGARv6_CH4_MonthlyScaleFactors_TNR_Aviation_CRS.0.1x0.1.nc
EDGAR_SEASONAL_SF_TNR_AV_LTO 1 2018 Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2023-04/EDGARv6_SF/EDGARv6_CH4_MonthlyScaleFactors_TNR_Aviation_LTO.0.1x0.1.nc
EDGAR_SEASONAL_SF_TNR_Other  1 2018 Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2023-04/EDGARv6_SF/EDGARv6_CH4_MonthlyScaleFactors_TNR_Other.0.1x0.1.nc
EDGAR_SEASONAL_SF_TNR_SHIP   1 2018 Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2023-04/EDGARv6_SF/EDGARv6_CH4_MonthlyScaleFactors_TNR_Ship.0.1x0.1.nc
EDGAR_SEASONAL_SF_TRO_noRes  1 2018 Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2023-04/EDGARv6_SF/EDGARv6_CH4_MonthlyScaleFactors_TRO_noRES.0.1x0.1.nc
EDGAR_SEASONAL_SF_WWT        1 2018 Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2023-04/EDGARv6_SF/EDGARv6_CH4_MonthlyScaleFactors_WWT.0.1x0.1.nc

#==============================================================================
# %%%%% Scale factors for CO emissions %%%%%
#
# If USE_CO_DATA is set to "false" in HEMCO_Config.rc,
# then you must comment out all of these entries below.
#==============================================================================

# --- HTAP annual scale factors ---
# Add entries here

# --- CEDS vertical partitioning ---
ENERGY_LEVS   1 N Y - none none g_energy   ./HcoDir/VerticalScaleFactors/v2021-05/gc_layers.nc
INDUSTRY_LEVS 1 N Y - none none g_industry ./HcoDir/VerticalScaleFactors/v2021-05/gc_layers.nc
SHIP_LEVS     1 N Y - none none cmv_c3     ./HcoDir/VerticalScaleFactors/v2021-05/gc_layers.nc

#==============================================================================
# %%%%% Scale factors for CO2 emissions %%%%%
#
# If USE_CO2_DATA is set to "false" in HEMCO_Config.rc,
# then you must comment out all of these entries below.
#==============================================================================

# --- National fossil fuel CO2 scale factors (Nassar et al, 2013) ---
CO2_DIURNAL 1 N    Y F%y4-%m2-%d2T%h2:00:00 none none diurnal_scale_factors ./HcoDir/CO2/v2015-04/FOSSIL/TIMES_diurnal_scale_factors.nc
CO2_WEEKLY  1 2006 Y F%y4-%m2-%d2T00:00:00  none none weekly_scale_factors  ./HcoDir/CO2/v2015-04/FOSSIL/TIMES_weekly_scale_factors.nc

# --- Domestic aviation surface correction factor ---
AVIATION_SURF_CORR      1       2004 Y F%y2-01-01T00:00:00  none none CO2 ./HcoDir/CO2/v2022-11/FOSSIL/Aviation_SurfCorr_SclFac.1x1.nc
OCEANCO2_SCALED_MONTHLY kg/m2/s N    Y F%y2-%m2-01T00:00:00 none none CO2 ./HcoDir/CO2/v2022-11/OCEAN/Scaled_Ocean_CO2_monthly.nc

###############################################################################
###
### Masks (update if HEMCO_Config.rc changes)
###
###############################################################################

# Use updated files (v2024-01); these are COARDS-compliant.
MEX_MASK          1 N Y - none none MASK ./HcoDir/MASKS/v2024-01/Mexico_Mask.001x001.nc
MEX_MASK_MIRROR   1 N Y - none none MASK ./HcoDir/MASKS/v2024-01/Mexico_Mask_Mirror.001x001.nc
CAN_MASK          1 N Y - none none MASK ./HcoDir/MASKS/v2024-01/Canada_Mask.001x001.nc
CAN_MASK_MIRROR   1 N Y - none none MASK ./HcoDir/MASKS/v2024-01/Canada_Mask_Mirror.001x001.nc
#
CHINA_MASK        1 N Y - none none MASK ./HcoDir/MASKS/v2018-09/China_mask.generic.1x1.nc
CONUS_MASK        1 N Y - none none MASK ./HcoDir/MASKS/v2018-09/CONUS_Mask.001x001.nc
CONUS_MASK_MIRROR 1 N Y - none none MASK ./HcoDir/MASKS/v2018-09/CONUS_Mask_Mirror.001x001.nc
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
