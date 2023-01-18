#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: HEMCO_Diagn.rc
#
# !DESCRIPTION: Configuration file for netCDF diagnostic output from HEMCO.
#\\
#\\
# !REMARKS:
#  Customized for the carbon simulation.
#EOP
#------------------------------------------------------------------------------
#BOC
# Name               Spec  ExtNr  Cat Hier Dim  OutUnit  LongName

###############################################################################
#####  CH4 sources                                                         #####
###############################################################################

EmisCH4_Total         CH4     -1   -1   -1   2  kg/m2/s  CH4_emissions_from_all_sectors
EmisCH4_Oil           CH4      0    1   -1   2  kg/m2/s  CH4_emissions_from_oil
EmisCH4_Gas           CH4      0    2   -1   2  kg/m2/s  CH4_emissions_from_gas
EmisCH4_Coal          CH4      0    3   -1   2  kg/m2/s  CH4_emissions_from_coal
EmisCH4_Livestock     CH4      0    4   -1   2  kg/m2/s  CH4_emissions_from_livestock
EmisCH4_Landfills     CH4      0    5   -1   2  kg/m2/s  CH4_emissions_from_landfills
EmisCH4_Wastewater    CH4      0    6   -1   2  kg/m2/s  CH4_emissions_from_wastewater
EmisCH4_Rice          CH4      0    7   -1   2  kg/m2/s  CH4_emissions_from_rice
EmisCH4_OtherAnth     CH4      0    8   -1   2  kg/m2/s  CH4_emissions_from_other_anthropogenic_sources
EmisCH4_BiomassBurn   CH4      0    9   -1   2  kg/m2/s  CH4_emissions_from_biomass_burning
EmisCH4_Wetlands      CH4      0   10   -1   2  kg/m2/s  CH4_emissions_from_wetlands
EmisCH4_Seeps         CH4      0   11   -1   2  kg/m2/s  CH4_emissions_from_geological_seeps
EmisCH4_Lakes         CH4      0   12   -1   2  kg/m2/s  CH4_emissions_from_lakes
EmisCH4_Termites      CH4      0   13   -1   2  kg/m2/s  CH4_emissions_from_termites
EmisCH4_SoilAbsorb    CH4      0   14   -1   2  kg/m2/s  CH4_emissions_from_soil_absorption

###############################################################################
#####  CO sources                                                         #####
###############################################################################
EmisCO_Total          CO      -1   -1   -1   3  kg/m2/s  CO_emission_flux_from_all_sectors
EmisCO_Aircraft       CO       0   20   -1   3  kg/m2/s  CO_emission_flux_from_aircraft
EmisCO_Anthro         CO       0    1   -1   3  kg/m2/s  CO_emission_flux_from_anthropogenic
EmisCO_BioBurn        CO       0    5   -1   2  kg/m2/s  CO_emission_flux_from_biomass_burning
EmisCO_Ship           CO       0   10   -1   2  kg/m2/s  CO_emission_flux_from_ships

###############################################################################
#####  CO2 sources                                                        #####
###############################################################################
EmisCO2_Total         CO2      0   -1   -1   3  kg/m2/s  CO2_total_emissions
EmisCO2_FossilFuel    CO2      0    1   -1   2  kg/m2/s  CO2_anthropogenic_emissions
EmisCO2_Ocean         CO2      0    2   -1   2  kg/m2/s  CO2_ocean_emissions
EmisCO2_BalBiosph     CO2      0    3   -1   2  kg/m2/s  CO2_balanced_biosphere
EmisCO2_NetTerrExch   CO2      0    5   -1   2  kg/m2/s  CO2_net_terrestrial_exchange
EmisCO2_Ship          CO2      0    6   -1   2  kg/m2/s  CO2_ship_emissions
EmisCO2_Aviation      CO2      0    7   -1   3  kg/m2/s  CO2_aviation_emissions
EmisCO2_CO2SurfCorr   CO2      0    8   -1   2  kg/m2/s  CO2_surface_correction_for_CO_oxidation

#NOTE: Biomass burning doesn't seem to work.  Keep here for now.
# Uncomment this if GFED is used
#EmisCO2_Biomass_GFED  CO2    111    1    1   2  kg/m2/s  CO2_biomass_burning_emissions_from_GFED4

# Uncomment this if FINN is used
#EmisCO2_Biomass_FINN  CO2    114    1    1   2  kg/m2/s  CO2_biomass_burning_emissions_from_FINN

###############################################################################
#####  OCS sources                                                        #####
###############################################################################
EmisOCS_Total         OCS      0   -1   -1   3  kg/m2/s  OCS_emission_flux_from_all_sectors
EmisOCS_Anthro        OCS      0    1   -1   2  kg/m2/s  OCS_emission_flux_from_anthropogenic
EmisOCS_Bioburn       OCS      0    2   -1   2  kg/m2/s  OCS_emission_flux_from_biomass_burning
EmisOCS_MissingOcean  OCS      0    3   -1   2  kg/m2/s  OCS_emission_flux_from_missing_ocean
EmisOCS_Ocean         OCS      0    4   -1   2  kg/m2/s  OCS_emission_flux_from_ocean

# Tagged CO species, uncomment if needed
#
################################################################################
######  COus sources                                                       #####
################################################################################
#EmisCOus_Total        COus    -1   -1   -1   3  kg/m2/s  COus_emission_flux_from_all_sectors
#EmisCOus_Aircraft     COus     0   20   -1   3  kg/m2/s  COus_emission_flux_from_aircraft
#EmisCOus_Anthro       COus     0    1   -1   3  kg/m2/s  COus_emission_flux_from_anthropogenic
#EmisCOus_BioBurn      COus     0    5   -1   2  kg/m2/s  COus_emission_flux_from_biomass_burning
#EmisCOus_Ship         COus     0   10   -1   2  kg/m2/s  COus_emission_flux_from_ships
#
################################################################################
######  COeur sources                                                      #####
################################################################################
#EmisCOeur_Total       COeur   -1   -1   -1   3  kg/m2/s  COeur_emission_flux_from_all_sectors
#EmisCOeur_Aircraft    COeur    0   20   -1   3  kg/m2/s  COeur_emission_flux_from_aircraft
#EmisCOeur_Anthro      COeur    0    1   -1   3  kg/m2/s  COeur_emission_flux_from_anthropogenic
#EmisCOeur_BioBurn     COeur    0    5   -1   2  kg/m2/s  COeur_emission_flux_from_biomass_burning
#EmisCOeur_Ship        COeur    0   10   -1   2  kg/m2/s  COeur_emission_flux_from_ships
#
################################################################################
######  COasia sources                                                     #####
################################################################################
#EmisCOasia_Total      COasia  -1   -1   -1   3  kg/m2/s  COasia_emission_flux_from_all_sectors
#EmisCOasia_Aircraft   COasia   0   20   -1   3  kg/m2/s  COasia_emission_flux_from_aircraft
#EmisCOasia_Anthro     COasia   0    1   -1   3  kg/m2/s  COasia_emission_flux_from_anthropogenic
#EmisCOasia_BioBurn    COasia   0    5   -1   2  kg/m2/s  COasia_emission_flux_from_biomass_burning
#EmisCOasia_Ship       COasia   0   10   -1   2  kg/m2/s  COasia_emission_flux_from_ships
#
################################################################################
######  COoth sources                                                      #####
################################################################################
#EmisCOoth_Total       COoth   -1   -1   -1   3  kg/m2/s  COoth_emission_flux_from_all_sectors
#EmisCOoth_Aircraft    COoth    0   20   -1   3  kg/m2/s  COoth_emission_flux_from_aircraft
#EmisCOoth_Anthro      COoth    0    1   -1   3  kg/m2/s  COoth_emission_flux_from_anthropogenic
#EmisCOoth_BioBurn     COoth    0    5   -1   2  kg/m2/s  COoth_emission_flux_from_biomass_burning
#EmisCOoth_Ship        COoth    0   10   -1   2  kg/m2/s  COoth_emission_flux_from_ships
#
################################################################################
######  COus sources                                                       #####
################################################################################
#EmisCOus_Total        COus    -1   -1   -1   3  kg/m2/s  COus_emission_flux_from_all_sectors
#EmisCOus_Aircraft     COus     0   20   -1   3  kg/m2/s  COus_emission_flux_from_aircraft
#EmisCOus_Anthro       COus     0    1   -1   3  kg/m2/s  COus_emission_flux_from_anthropogenic
#EmisCOus_BioBurn      COus     0    5   -1   2  kg/m2/s  COus_emission_flux_from_biomass_burning
#EmisCOus_Ship         COus     0   10   -1   2  kg/m2/s  COus_emission_flux_from_ships
#
################################################################################
######  COeur sources                                                      #####
################################################################################
#EmisCOeur_Total       COeur   -1   -1   -1   3  kg/m2/s  COeur_emission_flux_from_all_sectors
#EmisCOeur_Aircraft    COeur    0   20   -1   3  kg/m2/s  COeur_emission_flux_from_aircraft
#EmisCOeur_Anthro      COeur    0    1   -1   3  kg/m2/s  COeur_emission_flux_from_anthropogenic
#EmisCOeur_BioBurn     COeur    0    5   -1   2  kg/m2/s  COeur_emission_flux_from_biomass_burning
#EmisCOeur_Ship        COeur    0   10   -1   2  kg/m2/s  COeur_emission_flux_from_ships
#
################################################################################
######  COasia sources                                                     #####
################################################################################
#EmisCOasia_Total      COasia  -1   -1   -1   3  kg/m2/s  COasia_emission_flux_from_all_sectors
#EmisCOasia_Aircraft   COasia   0   20   -1   3  kg/m2/s  COasia_emission_flux_from_aircraft
#EmisCOasia_Anthro     COasia   0    1   -1   3  kg/m2/s  COasia_emission_flux_from_anthropogenic
#EmisCOasia_BioBurn    COasia   0    5   -1   2  kg/m2/s  COasia_emission_flux_from_biomass_burning
#EmisCOasia_Ship       COasia   0   10   -1   2  kg/m2/s  COasia_emission_flux_from_ships
#
################################################################################
######  COoth sources                                                      #####
################################################################################
#EmisCOoth_Total       COoth   -1   -1   -1   3  kg/m2/s  COoth_emission_flux_from_all_sectors
#EmisCOoth_Aircraft    COoth    0   20   -1   3  kg/m2/s  COoth_emission_flux_from_aircraft
#EmisCOoth_Anthro      COoth    0    1   -1   3  kg/m2/s  COoth_emission_flux_from_anthropogenic
#EmisCOoth_BioBurn     COoth    0    5   -1   2  kg/m2/s  COoth_emission_flux_from_biomass_burning
#EmisCOoth_Ship        COoth    0   10   -1   2  kg/m2/s  COoth_emission_flux_from_ships