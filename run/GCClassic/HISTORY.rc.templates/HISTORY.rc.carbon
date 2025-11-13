###############################################################################
###  HISTORY.rc file for GEOS-Chem aerosol-only specialty simulations       ###
###  Contact: GEOS-Chem Support Team (geos-chem-support@g.harvard.edu)      ###
###############################################################################

#============================================================================
# EXPID allows you to specify the beginning of the file path corresponding
# to each diagnostic collection.  For example:
#
#   EXPID: ./GEOSChem
#      Will create netCDF files whose names begin "GEOSChem",
#      in this run directory.
#
#   EXPID: ./OutputDir/GEOSChem
#      Will create netCDF files whose names begin with "GEOSChem"
#      in the OutputDir sub-folder of this run directory.
#
#============================================================================
EXPID:  ./OutputDir/GEOSChem

#==============================================================================
# %%%%% COLLECTION NAME DECLARATIONS %%%%%
#
# To enable a collection, remove the "#" character in front of its name. The
# Restart collection should always be turned on.
#
# NOTE: These are the "default" collections for GEOS-Chem, but you can create
# your own customized diagnostic collections as well.
#==============================================================================
COLLECTIONS: 'Restart',
             'SpeciesConc',
             #'Budget',
             #'CloudConvFlux',
             #'DryDep',
             'ProdLoss',
             #'SatDiagn',
             #'SatDiagnEdge',
             #'StateChm',
             'StateMet',
             #'StateMetLevEdge',
             #'BoundaryConditions',
::
###############################################################################
### The rest of this file consists of collection definitions.               ###
### Above collections whose declarations are commented out will be ignored. ###
### Make sure that each collection definition, as well as the list of       ###
### collections above, ends with a double-colon ("::").                     ###
###############################################################################

#==============================================================================
# %%%%% THE Restart COLLECTION %%%%%
#
# GEOS-Chem restart file fields
#==============================================================================
  Restart.filename:           './Restarts/GEOSChem.Restart.%y4%m2%d2_%h2%n2z.nc4',
  Restart.frequency:          'End',
  Restart.duration:           'End',
  Restart.mode:               'instantaneous'
  Restart.fields:             'SpeciesRst_?ALL?               ',
                              'Met_DELPDRY                    ',
                              'Met_BXHEIGHT                   ',
                              'Met_TropLev                    ',
::
#==============================================================================
# %%%%% THE SpeciesConc COLLECTION %%%%%
#
# GEOS-Chem species concentrations (default = all species)
#
# Concentrations may be saved out as SpeciesConcVV  [v/v dry air] or
#                                    SpeciesConcMND [molec/cm3]
#==============================================================================
  SpeciesConc.template:       '%y4%m2%d2_%h2%n2z.nc4',
  SpeciesConc.frequency:      ${RUNDIR_HIST_TIME_AVG_FREQ}
  SpeciesConc.duration:       ${RUNDIR_HIST_TIME_AVG_DUR}
  SpeciesConc.mode:           'time-averaged'
  SpeciesConc.fields:         'SpeciesConcVV_?ALL?           ',
                              #'SpeciesConcMND_?ALL?          ',
::
#==============================================================================
# %%%%% THE Budget COLLECTION %%%%%
#
# GEOS-Chem budget diagnostics defined as species kg/s in the column
# (full, troposphere, or PBL) due to a single component (e.g. chemistry)
# (default = advected species)
#==============================================================================
  Budget.template:     '%y4%m2%d2_%h2%n2z.nc4',
  Budget.frequency:    ${RUNDIR_HIST_TIME_AVG_FREQ}
  Budget.duration:     ${RUNDIR_HIST_TIME_AVG_DUR}
  Budget.mode:         'time-averaged'
  Budget.fields:       'BudgetEmisDryDepFull_?ADV?           ',
                       'BudgetEmisDryDepTrop_?ADV?           ',
                       'BudgetEmisDryDepPBL_?ADV?            ',
                       'BudgetEmisDryDepLevs1to35_?ADV?      ',
                       'BudgetChemistryFull_?ADV?            ',
                       'BudgetChemistryTrop_?ADV?            ',
                       'BudgetChemistryPBL_?ADV?             ',
                       'BudgetChemistryLevs1to35_?ADV?       ',
                       'BudgetTransportFull_?ADV?            ',
                       'BudgetTransportTrop_?ADV?            ',
                       'BudgetTransportPBL_?ADV?             ',
                       'BudgetTransportLevs1to35_?ADV?       ',
                       'BudgetMixingFull_?ADV?               ',
                       'BudgetMixingTrop_?ADV?               ',
                       'BudgetMixingPBL_?ADV?                ',
                       'BudgetMixingLevs1to35_?ADV?          ',
                       'BudgetConvectionFull_?ADV?           ',
                       'BudgetConvectionTrop_?ADV?           ',
                       'BudgetConvectionPBL_?ADV?            ',
                       'BudgetConvectionLevs1to35_?ADV?      ',
                       'BudgetWetDepFull_?WET?               ',
                       'BudgetWetDepTrop_?WET?               ',
                       'BudgetWetDepPBL_?WET?                ',
                       'BudgetWetDepLevs1to35_?WET?          ',
::
#==============================================================================
# %%%%% THE CloudConvFlux COLLECTION %%%%%
#
# Cloud convective flux (default = advected species)
#==============================================================================
  CloudConvFlux.template:     '%y4%m2%d2_%h2%n2z.nc4',
  CloudConvFlux.frequency:    ${RUNDIR_HIST_TIME_AVG_FREQ}
  CloudConvFlux.duration:     ${RUNDIR_HIST_TIME_AVG_DUR}
  CloudConvFlux.mode:         'time-averaged'
  CloudConvFlux.fields:       'CloudConvFlux_?ADV?           ',
::
#==============================================================================
# %%%%% THE DryDep COLLECTION %%%%%
#
# Dry depositon fluxes and velocities (for all dry-depositing species)
#==============================================================================
  DryDep.template:            '%y4%m2%d2_%h2%n2z.nc4',
  DryDep.frequency:           ${RUNDIR_HIST_TIME_AVG_FREQ}
  DryDep.duration:            ${RUNDIR_HIST_TIME_AVG_DUR}
  DryDep.mode:                'time-averaged'
  DryDep.fields:              'DryDepVel_?DRY?               ',
                              'DryDep_?DRY?                  ',
                              #'DryDepChm_?DRY?              ',
                              #'DryDepMix_?DRY?              ',
::
#==============================================================================
# %%%%% THE ProdLoss COLLECTION %%%%%
#
# Production and loss fields from the carbon simulation via KPP
#==============================================================================
  ProdLoss.template:         '%y4%m2%d2_%h2%n2z.nc4',
  ProdLoss.frequency:        ${RUNDIR_HIST_TIME_AVG_FREQ}
  ProdLoss.duration:         ${RUNDIR_HIST_TIME_AVG_DUR}
  ProdLoss.mode:             'time-averaged'
  ProdLoss.fields:           'ProdCOfromCH4                  ',
                             'ProdCOfromNMVOC                ',
                             'ProdCO2fromCO                  ',
::
#==============================================================================
# %%%%% THE SatDiagn COLLECTION %%%%%
#
# GEOS-Chem data during satellite overpass
#==============================================================================
  SatDiagn.template:          '%y4%m2%d2_%h2%n2z.nc4',
  SatDiagn.format:            'CFIO',
  SatDiagn.frequency:         00000001 000000
  SatDiagn.duration:          00000100 000000
  SatDiagn.hrrange:           11.98 15.02
  SatDiagn.mode:              'time-averaged'
  SatDiagn.fields:            'SatDiagnConc_BCPI              ',
                              'SatDiagnConc_BCPO              ',
                              'SatDiagnRH                     ',
                              'SatDiagnAirDen                 ',
                              'SatDiagnBoxHeight              ',
                              'SatDiagnPEdge                  ',
                              'SatDiagnTROPP                  ',
                              'SatDiagnPBLHeight              ',
                              'SatDiagnPBLTop                 ',
                              'SatDiagnPMID                   ',
                              'SatDiagnTAir                   ',
                              'SatDiagnSLP                    ',
                              'SatDiagnSPHU                   ',
                              'SatDiagnTS                     ',
                              'SatDiagnPBLTOPL                ',
                              'SatDiagnMODISLAI               ',
                              'SatDiagnWetLossLS_BCPI         ',
                              'SatDiagnWetLossConv_BCPI       ',
                              'SatDiagnDryDep_BCPI            ',
                              'SatDiagnDryDepVel_BCPI         ',
                              'SatDiagnColEmis_BCPI           ',
                              'SatDiagnSurfFlux_BCPI          ',
                              'SatDiagnProd_?PRD?             ',
                              'SatDiagnLoss_?LOS?             ',
::
#==============================================================================
# %%%%% THE SatDiagnEdge COLLECTION %%%%%
#
# GEOS-Chem data (on level edges) during satellite overpass
#==============================================================================
  SatDiagnEdge.template:      '%y4%m2%d2_%h2%n2z.nc4',
  SatDiagnEdge.format:        'CFIO',
  SatDiagnEdge.frequency:     00000001 000000
  SatDiagnEdge.duration:      00000100 000000
  SatDiagnEdge.hrrange:       11.98 15.02
  SatDiagnEdge.mode:          'time-averaged'
  SatDiagnEdge.fields:        'SatDiagnPEDGE                  ',
::
#==============================================================================
# %%%%% THE StateChm COLLECTION %%%%%
#
# Fields of the State_Chm object (see also Aerosols collection)
#
# Available for full-chemistry and aerosols-only simulations
#==============================================================================
  StateChm.template:          '%y4%m2%d2_%h2%n2z.nc4',
  StateChm.frequency:         ${RUNDIR_HIST_TIME_AVG_FREQ}
  StateChm.duration:          ${RUNDIR_HIST_TIME_AVG_DUR}
  StateChm.mode:              'time-averaged'
  StateChm.fields:            'Chem_IsorropAeropHAccum       ',
                              'Chem_IsorropAeropHCoarse      ',
                              'Chem_IsorropHplusAccum        ',
                              'Chem_IsorropHplusCoarse       ',
                              'Chem_IsorropAeroH2OAccum      ',
                              'Chem_IsorropAeroH2OCoarse     ',
                              'Chem_IsorropSulfate           ',
                              'Chem_IsorropNitrateAccum      ',
                              'Chem_IsorropNitrateCoarse     ',
                              'Chem_IsorropChlorideAccum     ',
                              'Chem_IsorropChlorideCoarse    ',
                              'Chem_IsorropBisulfate         ',
                              'Chem_pHCloud                  ',
                              'Chem_isCloud                  ',
                              'Chem_SSAlkAccumMode           ',
                              'Chem_SSAlkCoarseMode          ',
                              'Chem_HSO3AQ                   ',
                              'Chem_SO3AQ                    ',
                              'Chem_fupdateHOBr              ',
                              'Chem_GammaN2O5overall         ',
                              'Chem_GammaN2O5fine            ',
                              'Chem_YieldClNO2fine           ',
::
#==============================================================================
# %%%%% The StateMet COLLECTION %%%%%
#
# Fields of the State_Met object (also see the StateMetLevEdge collection)
#==============================================================================
  StateMet.template:          '%y4%m2%d2_%h2%n2z.nc4',
  StateMet.frequency:         ${RUNDIR_HIST_TIME_AVG_FREQ}
  StateMet.duration:          ${RUNDIR_HIST_TIME_AVG_DUR}
  StateMet.mode:              'time-averaged'
  StateMet.fields:            'Met_AD                        ',
                              'Met_AIRDEN                    ',
                              'Met_AIRVOL                    ',
                              'Met_ALBD                      ',
                              'Met_AVGW                      ',
                              'Met_BXHEIGHT                  ',
                              'Met_ChemGridLev               ',
                              'Met_CLDF                      ',
                              'Met_CLDFRC                    ',
                              'Met_CLDTOPS                   ',
                              'Met_DELP                      ',
                              'Met_DQRCU                     ',
                              'Met_DQRLSAN                   ',
                              'Met_DTRAIN                    ',
                              'Met_EFLUX                     ',
                              'Met_FRCLND                    ',
                              'Met_FRLAKE                    ',
                              'Met_FRLAND                    ',
                              'Met_FRLANDICE                 ',
                              'Met_FROCEAN                   ',
                              'Met_FRSEAICE                  ',
                              'Met_FRSNOW                    ',
                              'Met_GWETROOT                  ',
                              'Met_GWETTOP                   ',
                              'Met_HFLUX                     ',
                              'Met_LAI                       ',
                              'Met_PARDR                     ',
                              'Met_PARDF                     ',
                              'Met_PBLTOPL                   ',
                              'Met_PBLH                      ',
                              'Met_PHIS                      ',
                              'Met_PMID                      ',
                              'Met_PMIDDRY                   ',
                              'Met_PRECANV                   ',
                              'Met_PRECCON                   ',
                              'Met_PRECLSC                   ',
                              'Met_PRECTOT                   ',
                              'Met_PS1DRY                    ',
                              'Met_PS1WET                    ',
                              'Met_PS2DRY                    ',
                              'Met_PS2WET                    ',
                              'Met_PSC2WET                   ',
                              'Met_PSC2DRY                   ',
                              'Met_QI                        ',
                              'Met_QL                        ',
                              'Met_OMEGA                     ',
                              'Met_OPTD                      ',
                              'Met_REEVAPCN                  ',
                              'Met_REEVAPLS                  ',
                              'Met_SLP                       ',
                              'Met_SNODP                     ',
                              'Met_SNOMAS                    ',
                              'Met_SPHU                      ',
                              'Met_SPHU1                     ',
                              'Met_SPHU2                     ',
                              'Met_SUNCOS                    ',
                              'Met_SUNCOSmid                 ',
                              'Met_SWGDN                     ',
                              'Met_T                         ',
                              'Met_TAUCLI                    ',
                              'Met_TAUCLW                    ',
                              'Met_THETA                     ',
                              'Met_TMPU1                     ',
                              'Met_TMPU2                     ',
                              'Met_TO3                       ',
                              'Met_TropHt                    ',
                              'Met_TropLev                   ',
                              'Met_TropP                     ',
                              'Met_TS                        ',
                              'Met_TSKIN                     ',
                              'Met_TV                        ',
                              'Met_U                         ',
                              'Met_U10M                      ',
                              'Met_USTAR                     ',
                              'Met_UVALBEDO                  ',
                              'Met_V                         ',
                              'Met_V10M                      ',
                              'Met_Z0                        ',
::
#==============================================================================
# %%%%% THE StateMetLevEdge COLLECTION %%%%%
#
# Diagnostics that are defined on grid box level edges
#==============================================================================
  StateMetLevEdge.template:    '%y4%m2%d2_%h2%n2z.nc4',
  StateMetLevEdge.frequency:   ${RUNDIR_HIST_TIME_AVG_FREQ}
  StateMetLevEdge.duration:    ${RUNDIR_HIST_TIME_AVG_DUR}
  StateMetLevEdge.mode:        'time-averaged'
  StateMetLevEdge.fields:      'Met_CMFMC                     ',
                               'Met_PEDGE                     ',
                               'Met_PEDGEDRY                  ',
                               'Met_PFICU                     ',
                               'Met_PFILSAN                   ',
                               'Met_PFLCU                     ',
                               'Met_PFLLSAN                   ',
::
#==============================================================================
# %%%%% THE BoundaryConditions COLLECTION %%%%%
#
# GEOS-Chem boundary conditions for use in nested grid simulations
#
# Available for all simulations
#
# NOTE: Do not change the BoundaryConditions.frequency setting below,
# because GEOS-Chem nested-grid simulations expect to read boundary
# condition data at 3 hour intervals.
#==============================================================================
  BoundaryConditions.template:   '%y4%m2%d2_%h2%n2z.nc4',
  BoundaryConditions.frequency:  00000000 030000
  BoundaryConditions.duration:   00000001 000000
  BoundaryConditions.mode:       'instantaneous'
  BoundaryConditions.fields:     'SpeciesBC_?ADV?             ',
::