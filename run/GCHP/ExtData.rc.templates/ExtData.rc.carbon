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

#==============================================================================
# --- Fields for lightning emissions (Extension 103) ---
# Not used but included in Chem_Registry.rc import list
#==============================================================================
FLASH_DENS 1 N N - none none LDENS /dev/null
CONV_DEPTH 1 N N - none none CTH   /dev/null

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

# --- Gridded GHGI v2 Express Extension (Maasakkers et al., ES&T, 2023) ---
# Include scaling to convert molecules/cm2/s to kg/m2/s
GHGI_EE_OIL_EXPLORATION   kg/m2/s N Y F%y4-01-01T00:00:00 none 2.66350462E-22 emi_ch4_1B2a_Petroleum_Systems_Exploration   ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_OIL_PRODUCTION    kg/m2/s N Y F%y4-01-01T00:00:00 none 2.66350462E-22 emi_ch4_1B2a_Petroleum_Systems_Production    ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_OIL_REFINING      kg/m2/s N Y F%y4-01-01T00:00:00 none 2.66350462E-22 emi_ch4_1B2a_Petroleum_Systems_Refining      ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_OIL_TRANSPORT     kg/m2/s N Y F%y4-01-01T00:00:00 none 2.66350462E-22 emi_ch4_1B2a_Petroleum_Systems_Transport     ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_GAS_DISTRIBUTION  kg/m2/s N Y F%y4-01-01T00:00:00 none 2.66350462E-22 emi_ch4_1B2b_Natural_Gas_Distribution        ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_GAS_EXPLORATION   kg/m2/s N Y F%y4-01-01T00:00:00 none 2.66350462E-22 emi_ch4_1B2b_Natural_Gas_Exploration         ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_GAS_PROCESSING    kg/m2/s N Y F%y4-01-01T00:00:00 none 2.66350462E-22 emi_ch4_1B2b_Natural_Gas_Processing          ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_GAS_PRODUCTION    kg/m2/s N Y F%y4-01-01T00:00:00 none 2.66350462E-22 emi_ch4_1B2b_Natural_Gas_Production          ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_GAS_TRANSMISSION  kg/m2/s N Y F%y4-01-01T00:00:00 none 2.66350462E-22 emi_ch4_1B2b_Natural_Gas_TransmissionStorage ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_GAS_POSTMETER     kg/m2/s N Y F%y4-01-01T00:00:00 none 2.66350462E-22 emi_ch4_Supp_1B2b_PostMeter                  ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAL_UNDERGROUND  kg/m2/s N Y F%y4-01-01T00:00:00 none 2.66350462E-22 emi_ch4_1B1a_Underground_Coal                ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAL_SURFACE      kg/m2/s N Y F%y4-01-01T00:00:00 none 2.66350462E-22 emi_ch4_1B1a_Surface_Coal                    ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_COAL_ABANDONED    kg/m2/s N Y F%y4-01-01T00:00:00 none 2.66350462E-22 emi_ch4_1B1a_Abandoned_Coal                  ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_LIVESTOCK_ENT     kg/m2/s N Y F%y4-01-01T00:00:00 none 2.66350462E-22 emi_ch4_3A_Enteric_Fermentation              ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_LIVESTOCK_MAN     kg/m2/s N Y F%y4-01-01T00:00:00 none 2.66350462E-22 emi_ch4_3B_Manure_Management                 ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_LANDFILLS_IND     kg/m2/s N Y F%y4-01-01T00:00:00 none 2.66350462E-22 emi_ch4_5A1_Landfills_Industrial             ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_LANDFILLS_MSW     kg/m2/s N Y F%y4-01-01T00:00:00 none 2.66350462E-22 emi_ch4_5A1_Landfills_MSW                    ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_LANDFILLS_COMP    kg/m2/s N Y F%y4-01-01T00:00:00 none 2.66350462E-22 emi_ch4_5B1_Composting                       ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_WASTEWATER_DOM    kg/m2/s N Y F%y4-01-01T00:00:00 none 2.66350462E-22 emi_ch4_5D_Wastewater_Treatment_Domestic     ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_WASTEWATER_IND    kg/m2/s N Y F%y4-01-01T00:00:00 none 2.66350462E-22 emi_ch4_5D_Wastewater_Treatment_Industrial   ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_RICE              kg/m2/s N Y F%y4-01-01T00:00:00 none 2.66350462E-22 emi_ch4_3C_Rice_Cultivation                  ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_OTHER_MCOMB       kg/m2/s N Y F%y4-01-01T00:00:00 none 2.66350462E-22 emi_ch4_1A_Combustion_Mobile                 ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_OTHER_SCOMB       kg/m2/s N Y F%y4-01-01T00:00:00 none 2.66350462E-22 emi_ch4_1A_Combustion_Stationary             ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_OTHER_PIND        kg/m2/s N Y F%y4-01-01T00:00:00 none 2.66350462E-22 emi_ch4_2B8_Industry_Petrochemical           ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_OTHER_FIND        kg/m2/s N Y F%y4-01-01T00:00:00 none 2.66350462E-22 emi_ch4_2C2_Industry_Ferroalloy              ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_OTHER_BURN        kg/m2/s N Y F%y4-01-01T00:00:00 none 2.66350462E-22 emi_ch4_3F_Field_Burning                     ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc
GHGI_EE_OTHER_ABOG        kg/m2/s N Y F%y4-01-01T00:00:00 none 2.66350462E-22 emi_ch4_1B2ab_Abandoned_Oil_Gas              ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Express_Extension_Gridded_GHGI_Methane_v2_%y4.nc

# --- Scarpelli Mexico ---
# Include scaling to convert molecules/cm2/s to kg/m2/s
MEX_OIL          kg/m2/s N Y - none 2.66350462E-22 emis_ch4  ./HcoDir/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_oil_2015.nc
MEX_GAS          kg/m2/s N Y - none 2.66350462E-22 emis_ch4  ./HcoDir/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_gas_2015.nc
MEX_COAL         kg/m2/s N Y - none 2.66350462E-22 emis_ch4  ./HcoDir/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_coal_2015.nc
MEX_LIVESTOCK_A  kg/m2/s N Y - none 2.66350462E-22 emis_ch4  ./HcoDir/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_livestock_A_2015.nc
MEX_LIVESTOCK_B  kg/m2/s N Y - none 2.66350462E-22 emis_ch4  ./HcoDir/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_livestock_B_2015.nc
MEX_LANDFILLS    kg/m2/s N Y - none 2.66350462E-22 emis_ch4  ./HcoDir/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_landfill_2015.nc
MEX_WASTEWATER   kg/m2/s N Y - none 2.66350462E-22 emis_ch4  ./HcoDir/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_waste_2015.nc
MEX_RICE         kg/m2/s N Y - none 2.66350462E-22 emis_ch4  ./HcoDir/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_rice_2015.nc
MEX_OTHER        kg/m2/s N Y - none 2.66350462E-22 emis_ch4  ./HcoDir/CH4/v2022-11/Scarpelli_Mexico/MEX_Tia2020_other_anthro_2015.nc

# --- Scarpelli Canada ---
CAN_OIL_GAS_COMBUSTION  kg/m2/s N Y - none none  oil_gas_combustion_total   ./HcoDir/CH4/v2022-01/Scarpelli_Canada/can_emis_oil_gas_combustion_2018.nc
CAN_OIL_GAS_LEAKAGE     kg/m2/s N Y - none none  oil_gas_leakage_total      ./HcoDir/CH4/v2022-01/Scarpelli_Canada/can_emis_oil_gas_leakage_2018.nc
CAN_OIL_GAS_VENT_FLARE  kg/m2/s N Y - none none  oil_gas_vent_flare_total   ./HcoDir/CH4/v2022-01/Scarpelli_Canada/can_emis_oil_gas_vent_flare_2018.nc
CAN_COAL                kg/m2/s N Y - none none  coal_total                 ./HcoDir/CH4/v2022-01/Scarpelli_Canada/can_emis_coal_2018.nc
CAN_LIVESTOCK           kg/m2/s N Y - none none  livestock_total            ./HcoDir/CH4/v2022-01/Scarpelli_Canada/can_emis_livestock_2018.nc
CAN_SOLID_WASTE         kg/m2/s N Y - none none  solid_waste_total          ./HcoDir/CH4/v2022-01/Scarpelli_Canada/can_emis_solid_waste_2018.nc
CAN_WASTEWATER          kg/m2/s N Y - none none  wastewater_total           ./HcoDir/CH4/v2022-01/Scarpelli_Canada/can_emis_wastewater_2018.nc
CAN_OTHER               kg/m2/s N Y - none none  other_minor_sources_total  ./HcoDir/CH4/v2022-01/Scarpelli_Canada/can_emis_other_minor_sources_2018.nc

# --- GFEI ---
# Include scaling to convert molecules/cm2/s to kg/m2/s
GFEI_CH4_OIL  kg/m2/s N Y - none 2.66350462E-22 emis_ch4  ./HcoDir/CH4/v2022-01/GFEIv2/Global_Fuel_Exploitation_Inventory_v2_2019_Oil_All.nc
GFEI_CH4_GAS  kg/m2/s N Y - none 2.66350462E-22 emis_ch4  ./HcoDir/CH4/v2022-01/GFEIv2/Global_Fuel_Exploitation_Inventory_v2_2019_Gas_All.nc
GFEI_CH4_COAL kg/m2/s N Y - none 2.66350462E-22 emis_ch4  ./HcoDir/CH4/v2022-01/GFEIv2/Global_Fuel_Exploitation_Inventory_v2_2019_Coal.nc

# --- GRPI ---
GRPI_CH4_RICE kg/m2/s N Y F2022-%m2-01T00:00:00  none none emi_ch4 ./HcoDir/CH4/v2025-01/GRPI/GRPI_01x01.nc

# --- EDGAR v8.0 ---
EDGAR8_CH4_PRO_OIL             kg/m2/s N Y F%y4-%m2-01T00:00:00  none none emi_ch4 ./HcoDir/CH4/v2024-02/EDGARv8/%y4/v8.0_FT2022_GHG_CH4_%y4_PRO_OIL_flx.nc         
EDGAR8_CH4_REF_TRF             kg/m2/s N Y F%y4-%m2-01T00:00:00  none none emi_ch4 ./HcoDir/CH4/v2024-02/EDGARv8/%y4/v8.0_FT2022_GHG_CH4_%y4_REF_TRF_flx.nc         
EDGAR8_CH4_PRO_GAS             kg/m2/s N Y F%y4-%m2-01T00:00:00  none none emi_ch4 ./HcoDir/CH4/v2024-02/EDGARv8/%y4/v8.0_FT2022_GHG_CH4_%y4_PRO_GAS_flx.nc         
EDGAR8_CH4_PRO_COAL            kg/m2/s N Y F%y4-%m2-01T00:00:00  none none emi_ch4 ./HcoDir/CH4/v2024-02/EDGARv8/%y4/v8.0_FT2022_GHG_CH4_%y4_PRO_COAL_flx.nc        
##EDGAR8_CH4_ENF                 kg/m2/s N Y F%y4-%m2-01T00:00:00  none none emi_ch4 ./HcoDir/CH4/v2024-02/EDGARv8/%y4/v8.0_FT2022_GHG_CH4_%y4_ENF_flx.nc             
##EDGAR8_CH4_MNM                 kg/m2/s N Y F%y4-%m2-01T00:00:00  none none emi_ch4 ./HcoDir/CH4/v2024-02/EDGARv8/%y4/v8.0_FT2022_GHG_CH4_%y4_MNM_flx.nc             
EDGAR8_CH4_SWD_LDF             kg/m2/s N Y F%y4-%m2-01T00:00:00  none none emi_ch4 ./HcoDir/CH4/v2024-02/EDGARv8/%y4/v8.0_FT2022_GHG_CH4_%y4_SWD_LDF_flx.nc         
EDGAR8_CH4_WWT                 kg/m2/s N Y F%y4-%m2-01T00:00:00  none none emi_ch4 ./HcoDir/CH4/v2024-02/EDGARv8/%y4/v8.0_FT2022_GHG_CH4_%y4_WWT_flx.nc             
EDGAR8_CH4_AGS                 kg/m2/s N Y F%y4-%m2-01T00:00:00  none none emi_ch4 ./HcoDir/CH4/v2024-02/EDGARv8/%y4/v8.0_FT2022_GHG_CH4_%y4_AGS_flx.nc             
EDGAR8_CH4_ENE                 kg/m2/s N Y F%y4-%m2-01T00:00:00  none none emi_ch4 ./HcoDir/CH4/v2024-02/EDGARv8/%y4/v8.0_FT2022_GHG_CH4_%y4_ENE_flx.nc             
EDGAR8_CH4_IND                 kg/m2/s N Y F%y4-%m2-01T00:00:00  none none emi_ch4 ./HcoDir/CH4/v2024-02/EDGARv8/%y4/v8.0_FT2022_GHG_CH4_%y4_IND_flx.nc             
EDGAR8_CH4_TNR_Aviation_CDS    kg/m2/s N Y F%y4-%m2-01T00:00:00  none none emi_ch4 ./HcoDir/CH4/v2024-02/EDGARv8/%y4/v8.0_FT2022_GHG_CH4_%y4_TNR_Aviation_CDS_flx.nc
EDGAR8_CH4_TNR_Aviation_CRS    kg/m2/s N Y F%y4-%m2-01T00:00:00  none none emi_ch4 ./HcoDir/CH4/v2024-02/EDGARv8/%y4/v8.0_FT2022_GHG_CH4_%y4_TNR_Aviation_CRS_flx.nc
EDGAR8_CH4_TNR_Aviation_LTO    kg/m2/s N Y F%y4-%m2-01T00:00:00  none none emi_ch4 ./HcoDir/CH4/v2024-02/EDGARv8/%y4/v8.0_FT2022_GHG_CH4_%y4_TNR_Aviation_LTO_flx.nc
EDGAR8_CH4_TRO                 kg/m2/s N Y F%y4-%m2-01T00:00:00  none none emi_ch4 ./HcoDir/CH4/v2024-02/EDGARv8/%y4/v8.0_FT2022_GHG_CH4_%y4_TRO_flx.nc             
EDGAR8_CH4_TNR_Other           kg/m2/s N Y F%y4-%m2-01T00:00:00  none none emi_ch4 ./HcoDir/CH4/v2024-02/EDGARv8/%y4/v8.0_FT2022_GHG_CH4_%y4_TNR_Other_flx.nc       
EDGAR8_CH4_TNR_Ship            kg/m2/s N Y F%y4-%m2-01T00:00:00  none none emi_ch4 ./HcoDir/CH4/v2024-02/EDGARv8/%y4/v8.0_FT2022_GHG_CH4_%y4_TNR_Ship_flx.nc        
EDGAR8_CH4_RCO                 kg/m2/s N Y F%y4-%m2-01T00:00:00  none none emi_ch4 ./HcoDir/CH4/v2024-02/EDGARv8/%y4/v8.0_FT2022_GHG_CH4_%y4_RCO_flx.nc             
EDGAR8_CH4_CHE                 kg/m2/s N Y F%y4-%m2-01T00:00:00  none none emi_ch4 ./HcoDir/CH4/v2024-02/EDGARv8/%y4/v8.0_FT2022_GHG_CH4_%y4_CHE_flx.nc             
EDGAR8_CH4_IRO                 kg/m2/s N Y F%y4-%m2-01T00:00:00  none none emi_ch4 ./HcoDir/CH4/v2024-02/EDGARv8/%y4/v8.0_FT2022_GHG_CH4_%y4_IRO_flx.nc             

# --- EDGAR v7 livestock ---
EDGAR7_CH4_LIVESTOCK__4A      kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2023-04/EDGARv7/%y4/v7.0_FT2021_CH4_%y4_ENF.0.1x0.1.nc
EDGAR7_CH4_LIVESTOCK__4B      kg/m2/s N Y F%y4-%m2-01T00:00:00  none none  emi_ch4  ./HcoDir/CH4/v2023-04/EDGARv7/%y4/v7.0_FT2021_CH4_%y4_MNM.0.1x0.1.nc

## Comment out to avoid double counting with GFED
##EDGAR8_CH4_AWB                 kg/m2/s N Y F%y4-%m2-01T00:00:00  none none emi_ch4 ./HcoDir/CH4/v2024-02/EDGARv8/%y4/v8.0_FT2022_GHG_CH4_%y4_AWB_flx.nc             
EDGAR8_CH4_SWD_INC             kg/m2/s N Y F%y4-%m2-01T00:00:00  none none emi_ch4 ./HcoDir/CH4/v2024-02/EDGARv8/%y4/v8.0_FT2022_GHG_CH4_%y4_SWD_INC_flx.nc         

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

# --- JPL WetCHARTs v1.0 (Bloom et al., https://doi.org/10.3334/ORNLDAAC/1502) ---
# Include scaling to convert molecules/cm2/s to kg/m2/s
JPLW_CH4  kg/m2/s N Y F%y4-%m2-01T00:00:00 none 2.66350462E-22 emi_ch4  ./HcoDir/CH4/v2024-01/JPL_WetCharts/HEensemble/JPL_WetCharts_2010-2019.Ensemble_Mean.0.5x0.5.nc

# --- Geological seeps ---
CH4_SEEPS kg/m2/s N Y - none none  emi_ch4  ./HcoDir/CH4/v2020-04/Seeps/Etiope_CH4GeologicalEmis_ScaledToHmiel.1x1.nc

# --- Emissions from reservoirs (Delwich et al 2022) ---
CH4_RES_DAM kg/m2/s N Y - none none CH4emis  ./HcoDir/CH4/v2024-01/ResME/ResME_Dam_Emissions.0.1x0.1.nc
CH4_RES_SFC kg/m2/s N Y - none none CH4emis  ./HcoDir/CH4/v2024-01/ResME/ResME_Surface_Emissions.0.1x0.1.nc

# --- Emissions from termites (Fung et al 1991) ---
CH4_TERMITES    kg/m2/s N Y - none none  CH4  ./HcoDir/CH4/v2022-11/4x5/termites.geos.4x5.nc

# --- Soil absorption fro MeMo model (Murguia-Flores et al. 2018, GMD) ---
#
# Note: Use climatology by default
CH4_SOILABSORB kg/m2/s  2009 Y F%y4-%m2-01T00:00:00 none none  CH4uptake  ./HcoDir/CH4/v2019-10/MeMo_SoilAbs/MeMo_CH4uptake_Climatology.nc
#
# Annual emissions are disabled by default
#CH4_SOILABSORB kg/m2/s N    Y F%y4-%m2-01T00:00:00 none none  CH4uptake  ./HcoDir/CH4/v2019-10/MeMo_SoilAbs/MeMo_CH4uptake_1990-2009.nc

#==============================================================================
# %%%%% CO emissions %%%%%
#==============================================================================

# --- AEIC 2019 aircraft (AEIC) ---
AEIC19_MONMEAN_CO kg/m2/s Y Y F2019-%m2-01T00:00:00 none none CO ./HcoDir/AEIC2019/v2022-03/2019_monmean/AEIC_monmean_2019%m2.0.5x0.625.36L.nc
#
##AEIC19_DAILY_CO   kg/m2/s Y Y F%y4-%m2-%d2T00:00:00 none none CO ./HcoDir/AEIC2019/v2022-03/2019/%m2/AEIC_2019%m2%d2.0.5x0.625.36L.nc

# --- CEDSv2 emissions ---
CEDS_CO_AGR  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO_agr  ./HcoDir/CEDS/v2021-06/%y4/CO-em-anthro_CMIP_CEDS_%y4.nc
CEDS_CO_ENE  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO_ene  ./HcoDir/CEDS/v2021-06/%y4/CO-em-anthro_CMIP_CEDS_%y4.nc
CEDS_CO_IND  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO_ind  ./HcoDir/CEDS/v2021-06/%y4/CO-em-anthro_CMIP_CEDS_%y4.nc
CEDS_CO_TRA  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO_tra  ./HcoDir/CEDS/v2021-06/%y4/CO-em-anthro_CMIP_CEDS_%y4.nc
CEDS_CO_RCO  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO_rco  ./HcoDir/CEDS/v2021-06/%y4/CO-em-anthro_CMIP_CEDS_%y4.nc
CEDS_CO_SLV  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO_slv  ./HcoDir/CEDS/v2021-06/%y4/CO-em-anthro_CMIP_CEDS_%y4.nc
CEDS_CO_WST  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO_wst  ./HcoDir/CEDS/v2021-06/%y4/CO-em-anthro_CMIP_CEDS_%y4.nc
#
CEDS_CO_SHP  kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO_shp  ./HcoDir/CEDS/v2021-06/%y4/CO-em-anthro_CMIP_CEDS_%y4.nc

#==============================================================================
# --- CO2 emissions ---
#==============================================================================

# --- CO2 Fossil Fuel ---
FOSSILCO2_ODIAC  kg/m2/s N Y  F%y4-%m2-01T00:00:00 none none CO2 ./HcoDir/CO2/v2022-11/FOSSIL/ODIAC_CO2.monthly.generic.1x1.nc

# --- CO2 Ocean Exchange ---
OCEANCO2_SCALED_MONTHLY kg/m2/s N  Y F%y4-%m2-01T00:00:00 none none CO2 ./HcoDir/CO2/v2022-11/OCEAN/Scaled_Ocean_CO2_monthly.nc

# --- Balanced biosphere exchange ---
SIB_BBIO_CO2 kg/m2/s N Y F%y4-%m2-%d2T%h2:00:00 none none CO2 ./HcoDir/CO2/v2022-11/BIO/SiB3_3hr_NEP.nc

# --- Net Terrestrial Exchange ---
CO2_NET_TERRESTRIAL kg/m2/s N Y - none none CO2 ./HcoDir/CO2/v2022-11/BIO/Net_terrestrial_exch_5.29Pg.generic.1x1.nc

# --- CO2 Ship Emissions ---
CEDS_CO2_SHP kg/m2/s N Y F%y4-%m2-01T00:00:00 none none CO2_shp ./HcoDir/CEDS/v2021-06/%y4/CO2-em-anthro_CMIP_CEDS_%y4.nc

# --- CO2: AEIC 2019 aircraft emissions ---
AEIC19_MONMEAN_CO2 kg/m2/s Y Y F2019-%m2-01T00:00:00 none none FUELBURN ./HcoDir/AEIC2019/v2022-03/2019_monmean/AEIC_monmean_2019%m2.0.5x0.625.36L.nc

# --- CO2: Surface correction for CO oxidation ---
FOSSILCO2_MONTHLY kg/m/s N Y F%y4-%m2-01T01:00:00 none none CO2     ./HcoDir/CO2/v2022-11/FOSSIL/ODIAC_CO2.monthly.generic.1x1.nc
CO2_LIVESTOCK     kg/m/s Y Y F2004-%m2-01T01:00:00 none none CH4_004 ./HcoDir/CO2/v2022-11/CHEM/CH4_source.geos.2x25.nc
CO2_WASTE         kg/m/s Y Y F2004-%m2-01T01:00:00 none none CH4_005 ./HcoDir/CO2/v2022-11/CHEM/CH4_source.geos.2x25.nc
CO2_RICE          kg/m/s Y Y F2004-%m2-01T01:00:00 none none CH4_007 ./HcoDir/CO2/v2022-11/CHEM/CH4_source.geos.2x25.nc
CO2_WETLANDS      kg/m/s Y Y F2004-%m2-01T01:00:00 none none CH4_010 ./HcoDir/CO2/v2022-11/CHEM/CH4_source.geos.2x25.nc
CO2_NATURAL       kg/m/s Y Y F2004-%m2-01T01:00:00 none none CH4_012 ./HcoDir/CO2/v2022-11/CHEM/CH4_source.geos.2x25.nc
CO2_ISOPRENE      kg/m/s Y Y F2004-%m2-01T01:00:00 none none ISOP    ./HcoDir/CO2/v2022-11/CHEM/Isoprene-2004.geos.2x25.nc
CO2_MONOTERP      kg/m/s Y Y F2004-%m2-01T01:00:00 none none MONOT   ./HcoDir/CO2/v2022-11/CHEM/Monoterpene-2004.geos.2x25.nc

# --- CO2 production from CO ---
CO2_COPROD        kgC/m3/s N Y F%y4-%m2-01T00:00:00 none none LCO  ./HcoDir/CO2/v2024-01/CHEM/CO2_prod_rates.GEOS5.2x25.72L.nc

#==============================================================================
# --- OCS emission fluxes ---
#==============================================================================
OCS_ANTHRO        kgS/km2/s N Y F%y4-%m2-01T00:00:00 none 1.0e-6 COS_Flux ./HcoDir/OCS_FLUX/v2022-11/anthro_v3/v3_anthro/%y4/%m2.nc
OCS_BIOMASS       kgS/km2/s N Y F%y4-%m2-01T00:00:00 none 1.0e-6 COS_Flux ./HcoDir/OCS_FLUX/v2022-11/Biomassburn/CMSV2-CO2-GFED3-ZC/%y4/%m2.nc
OCS_MISSING_OCEAN kgS/km2/s N Y F%y4-%m2-01T00:00:00 none 1.0e-6 COS_Flux ./HcoDir/OCS_FLUX/v2022-11/MissingOcean/LUKAI/%y4/%m2.nc
OCS_OCEAN         kgS/km2/s N Y F%y4-%m2-01T00:00:00 none 1.0e-6 COS_Flux ./HcoDir/OCS_FLUX/v2022-11/OceanCOS/Kettle/%y4/%m2.nc

###############################################################################
###
### HEMCO Extension Data (update if HEMCO_Config.rc changes)
###
###############################################################################

# --- GFED4 (Extension 111)---
GFED_TEMP kgDM/m2/s N Y F%y4-%m2-01T00:00:00 none none DM_TEMP ./HcoDir/GFED4/v2023-03/%y4/GFED4_gen.025x025.%y4%m2.nc
GFED_AGRI kgDM/m2/s N Y F%y4-%m2-01T00:00:00 none none DM_AGRI ./HcoDir/GFED4/v2023-03/%y4/GFED4_gen.025x025.%y4%m2.nc
GFED_DEFO kgDM/m2/s N Y F%y4-%m2-01T00:00:00 none none DM_DEFO ./HcoDir/GFED4/v2023-03/%y4/GFED4_gen.025x025.%y4%m2.nc
GFED_BORF kgDM/m2/s N Y F%y4-%m2-01T00:00:00 none none DM_BORF ./HcoDir/GFED4/v2023-03/%y4/GFED4_gen.025x025.%y4%m2.nc
GFED_PEAT kgDM/m2/s N Y F%y4-%m2-01T00:00:00 none none DM_PEAT ./HcoDir/GFED4/v2023-03/%y4/GFED4_gen.025x025.%y4%m2.nc
GFED_SAVA kgDM/m2/s N Y F%y4-%m2-01T00:00:00 none none DM_SAVA ./HcoDir/GFED4/v2023-03/%y4/GFED4_gen.025x025.%y4%m2.nc

# --- GFED_daily scale factors ---
GFED_FRAC_DAY 1 N Y %y4-%m2-%d2T00:00:00 none none GFED_FRACDAY ./HcoDir/GFED4/v2023-03/%y4/GFED4_dailyfrac_gen.025x025.%y4%m2.nc

# --- GFED_3hourly --
#GFED_FRAC_3HOUR 1 N Y %y4-%m2-01T%h2:00:00 none none GFED_FRAC3HR ./HcoDir/GFED4/v2023-03/%y4/GFED4_3hrfrac_gen.025x025.%y4%m2.nc


###############################################################################
###
### HEMCO Non-Emissions Data (update if HEMCO_Config.rc changes)
###
###############################################################################

# --- Time zones (offset to UTC) ---
TIMEZONES count Y V 2017-%m2-01T00:00:00 none none UTC_OFFSET ./HcoDir/TIMEZONES/v2024-02/timezones_vohra_2017_0.1x0.1.nc

# --- OH from the latest 10-year benchmark ---
##GLOBAL_OH kg/m3 N N %y4-%m2-01T00:00:00       none none SpeciesConc_OH ./HcoDir/GCClassic_Output/14.0.0/%y4/GEOSChem.SpeciesConc.%y4%m201_0000z.nc4

# --- OH from GEOS-Chem v5-07-08 (used for CH4/IMI) ---
GLOBAL_OH kg/m3 1985 Y F%y4-%m2-01T00:00:00 none none OH ./HcoDir/OH/v2022-11/v5-07-08/OH_3Dglobal.geos5.72L.4x5.nc

# --- Global Cl concentrations ---
GLOBAL_Cl mol/mol N Y F%y4-%m2-01T00:00:00  none none SpeciesConc_Cl ./HcoDir/CH4/v2024-01/GCC_14_Output/%y4/GEOSChem.SpeciesConc.%y4%m201_0000z.nc4

# --- Stratospheric L(CO) from GMI ---
CH4_LOSS s-1 1985 Y F%y4-%m2-01T00:00:00 none none CH4loss ./HcoDir/CH4/v2024-01/GC_CH4_LOSS/GCC14_72LM.ch4loss.4x5.nc4

# --- P(CO) from CH4 and NMVOC from the last 10-yr benchmark ---
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
#==============================================================================

# --- GHGI seasonal scale factors ---
GHGI_OTH_STA_SF  1 N Y F%y4-%m2-01T00:00:00 none none monthly_scale_factor_1A_Combustion_Stationary           ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Gridded_GHGI_Methane_v2_Monthly_Scale_Factors_%y4.nc
GHGI_OIL_EXP_SF  1 N Y F%y4-%m2-01T00:00:00 none none monthly_scale_factor_1B2a_Petroleum_Systems_Exploration ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Gridded_GHGI_Methane_v2_Monthly_Scale_Factors_%y4.nc
GHGI_OIL_PRD_SF  1 N Y F%y4-%m2-01T00:00:00 none none monthly_scale_factor_1B2a_Petroleum_Systems_Production  ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Gridded_GHGI_Methane_v2_Monthly_Scale_Factors_%y4.nc
GHGI_OIL_REF_SF  1 N Y F%y4-%m2-01T00:00:00 none none monthly_scale_factor_1B2a_Petroleum_Systems_Refining    ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Gridded_GHGI_Methane_v2_Monthly_Scale_Factors_%y4.nc
GHGI_OIL_TRA_SF  1 N Y F%y4-%m2-01T00:00:00 none none monthly_scale_factor_1B2a_Petroleum_Systems_Transport   ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Gridded_GHGI_Methane_v2_Monthly_Scale_Factors_%y4.nc 
GHGI_GAS_EXP_SF  1 N Y F%y4-%m2-01T00:00:00 none none monthly_scale_factor_1B2b_Natural_Gas_Exploration       ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Gridded_GHGI_Methane_v2_Monthly_Scale_Factors_%y4.nc 
GHGI_GAS_PRD_SF  1 N Y F%y4-%m2-01T00:00:00 none none monthly_scale_factor_1B2b_Natural_Gas_Production        ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Gridded_GHGI_Methane_v2_Monthly_Scale_Factors_%y4.nc
GHGI_LIV_MAN_SF  1 N Y F%y4-%m2-01T00:00:00 none none monthly_scale_factor_3B_Manure_Management               ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Gridded_GHGI_Methane_v2_Monthly_Scale_Factors_%y4.nc 
GHGI_RIC_CUL_SF  1 N Y F%y4-%m2-01T00:00:00 none none monthly_scale_factor_3C_Rice_Cultivation                ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Gridded_GHGI_Methane_v2_Monthly_Scale_Factors_%y4.nc
GHGI_OTH_BUR_SF  1 N Y F%y4-%m2-01T00:00:00 none none monthly_scale_factor_3F_Field_Burning                   ./HcoDir/CH4/v2023-07/Gridded_GHGI_v2/Gridded_GHGI_Methane_v2_Monthly_Scale_Factors_%y4.nc

# --- Scarpelli_Mexico manure and rice scale factors ---
MANURE_SF  1 N    Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2017-10/Seasonal_SF/EMICH4_Manure_ScalingFactors.WithClimatology.nc
RICE_SF    1 2012 Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2017-10/Seasonal_SF/EMICH4_Rice_ScalingFactors.SetMissing.nc

# --- Seasonal scale factors for EDGARv7 ---
EDGAR_SEASONAL_SF_ENF 1 2018 Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2023-04/EDGARv6_SF/EDGARv6_CH4_MonthlyScaleFactors_ENF.0.1x0.1.nc
EDGAR_SEASONAL_SF_MNM 1 2018 Y F%y4-%m2-01T00:00:00 none none sf_ch4 ./HcoDir/CH4/v2023-04/EDGARv6_SF/EDGARv6_CH4_MonthlyScaleFactors_MNM.0.1x0.1.nc

#==============================================================================
# %%%%% Scale factors for CO emissions %%%%%
#==============================================================================

# --- CEDS vertical partitioning ---
ENERGY_LEVS   1 N Y - none none g_energy   ./HcoDir/VerticalScaleFactors/v2021-05/gc_layers.nc
INDUSTRY_LEVS 1 N Y - none none g_industry ./HcoDir/VerticalScaleFactors/v2021-05/gc_layers.nc
SHIP_LEVS     1 N Y - none none cmv_c3     ./HcoDir/VerticalScaleFactors/v2021-05/gc_layers.nc

#==============================================================================
# %%%%% Scale factors for CO2 emissions %%%%%
#==============================================================================

# --- National fossil fuel CO2 scale factors (Nassar et al, 2013) ---
CO2_DIURNAL 1 N Y F2006-01-01T%h2:00:00 none none diurnal_scale_factors ./HcoDir/CO2/v2015-04/FOSSIL/TIMES_diurnal_scale_factors.nc
CO2_WEEKLY  1 D Y F2006-01-%d2T00:00:00 none none weekly_scale_factors  ./HcoDir/CO2/v2015-04/FOSSIL/TIMES_weekly_scale_factors.nc

# --- Domestic aviation surface correction factor ---
AVIATION_SURF_CORR      1       Y  Y F2004-01-01T00:00:00  none none CO2 ./HcoDir/CO2/v2022-11/FOSSIL/Aviation_SurfCorr_SclFac.1x1.nc

###############################################################################
###
### Masks (update if HEMCO_Config.rc changes)
###
###############################################################################

# Comment these lines if you are not using Scarpelli_Mexico
MEX_MASK          1 N Y - none none MASK ./HcoDir/MASKS/v2024-04/Mexico_Mask.01x01.nc
MEX_MASK_MIRROR   1 N Y - none none MASK ./HcoDir/MASKS/v2024-04/Mexico_Mask_Mirror.01x01.nc

# Comment these lines if you are not using Scarpelli_Canada
CAN_MASK          1 N Y - none none MASK ./HcoDir/MASKS/v2024-04/Canada_Mask.01x01.nc
CAN_MASK_MIRROR   1 N Y - none none MASK ./HcoDir/MASKS/v2024-04/Canada_Mask_Mirror.01x01.nc

# Comment these lines if you are not using GHGI_v2/GHGIv2_Express_Ext
CONUS_MASK        1 N Y - none none MASK ./HcoDir/MASKS/v2024-04/CONUS_Mask.01x01.nc
CONUS_MASK_MIRROR 1 N Y - none none MASK ./HcoDir/MASKS/v2024-04/CONUS_Mask_Mirror.01x01.nc

# Comment these lines if you are not using any of the CEDS inventories
CHINA_MASK        1 N Y - none none MASK ./HcoDir/MASKS/v2018-09/China_mask.generic.1x1.nc

# For CH4 reservoir emissions from ResME
SEASONAL_RES_MASK 1 Y Y F2022-%m2-01T00:00:00 none none Mask ./HcoDir/CH4/v2024-01/ResME/reservoir_mask.01x01.nc
%%

DerivedExports%%
# ---------|---------|--------------------------------------------|
#  Export  | Primary |_________________ Mask _____________________|
#  Name    |  Name   |    Name    |           Expression          |
# ---------|---------|------------|-------------------------------|
# ---------|---------|------------|-------------------------------|
${RUNDIR_MET_EXTDATA_DERIVED_EXPORTS}
%%
