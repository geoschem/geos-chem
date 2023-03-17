###############################################################################
###  HISTORY.rc file for GEOS-Chem carbon simulation                   ###
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
             'Metrics',
             'SpeciesConc',
             #'Budget',
	     #'Carbon',
             #'CloudConvFlux',
             #'ConcAfterChem',
             #'LevelEdgeDiags',
             #'StateMet',
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
#
# Available for all simulations
#
# NOTE: Chem_KPPHvalue isn't needed since the carbon simulation uses
# the forward Euler solver instead of Rosenbrock.
#==============================================================================
  Restart.filename:     './Restarts/GEOSChem.Restart.%y4%m2%d2_%h2%n2z.nc4',
  Restart.frequency:    ${RUNDIR_HIST_INST_FREQ}
  Restart.duration:     ${RUNDIR_HIST_INST_DUR}
  Restart.mode:         'instantaneous'
  Restart.fields:       'SpeciesRst_?ALL?                     ',
                        'Met_DELPDRY                          ',
                        'Met_BXHEIGHT                         ',
                        'Met_TropLev                          ',
::
#==============================================================================
# %%%%% THE Metrics COLLECTION %%%%%
#
# Diagnostics for chemistry metrics such as global mean OH concentration,
# MCF lifetime, and CH4 lifetime.
#
# This diagnostic collection should always be left on.
# Use the run-directory script ./metrics_ch4.py to print results.
#
# Available for full-chemistry and CH4 simulations only.
#==============================================================================
  Metrics.template:           '%y4%m2%d2_%h2%n2z.nc4',
  Metrics.frequency:          'End',
  Metrics.duration:           'End',
  Metrics.mode:               'time-averaged'
  Metrics.fields:             'AirMassColumnFull             ',
                              'CH4emission                   ',
                              'CH4massColumnFull             ',
                              'CH4massColumnTrop             ',
                              'LossOHbyCH4columnTrop         ',
                              'LossOHbyMCFcolumnTrop         ',
                              'OHwgtByAirMassColumnFull      ',
::
#==============================================================================
# %%%%% THE SpeciesConc COLLECTION %%%%%
#
# GEOS-Chem species concentrations (default = advected species)
#
# Available for all simulations
#==============================================================================
  SpeciesConc.template:       '%y4%m2%d2_%h2%n2z.nc4',
  SpeciesConc.frequency:      ${RUNDIR_HIST_TIME_AVG_FREQ}
  SpeciesConc.duration:       ${RUNDIR_HIST_TIME_AVG_DUR}
  SpeciesConc.mode:           'time-averaged'
  SpeciesConc.fields:         'SpeciesConcVV_?ADV?           ',
                              'SpeciesConcMND_?ADV?          ',
::
#==============================================================================
# %%%%% THE Budget COLLECTION %%%%%
#
# GEOS-Chem budget diagnostics defined as species kg/s in the column
# (full, troposphere, or PBL) due to a single component (e.g. chemistry)
# (default = advected species)
#
# Available for all simulations
#==============================================================================
  Budget.template:            '%y4%m2%d2_%h2%n2z.nc4',
  Budget.frequency:           ${RUNDIR_HIST_TIME_AVG_FREQ}
  Budget.duration:            ${RUNDIR_HIST_TIME_AVG_DUR}
  Budget.mode:                'time-averaged'
  Budget.fields:              'BudgetEmisDryDepFull_?ADV?    ',
                              'BudgetEmisDryDepTrop_?ADV?    ',
                              'BudgetEmisDryDepPBL_?ADV?     ',
                              'BudgetChemistryFull_?ADV?     ',
                              'BudgetChemistryTrop_?ADV?     ',
                              'BudgetChemistryPBL_?ADV?      ',
                              'BudgetTransportFull_?ADV?     ',
                              'BudgetTransportTrop_?ADV?     ',
                              'BudgetTransportPBL_?ADV?      ',
                              'BudgetMixingFull_?ADV?        ',
                              'BudgetMixingTrop_?ADV?        ',
                              'BudgetMixingPBL_?ADV?         ',
                              'BudgetConvectionFull_?ADV?    ',
                              'BudgetConvectionTrop_?ADV?    ',
                              'BudgetConvectionPBL_?ADV?     ',
                              'BudgetWetDepFull_?WET?        ',
                              'BudgetWetDepTrop_?WET?        ',
                              'BudgetWetDepPBL_?WET?         ',
::
#==============================================================================
# %%%%% THE Carbon COLLECTION %%%%%
#
# Production and loss fields from the carbon simulation via KPP
#==============================================================================
  Carbon.template:       '%y4%m2%d2_%h2%n2z.nc4',
  Carbon.frequency:      ${RUNDIR_HIST_TIME_AVG_FREQ}
  Carbon.duration:       ${RUNDIR_HIST_TIME_AVG_DUR}
  Carbon.mode:           'time-averaged'
  Carbon.fields:         'OHconcAfterChem               ',
                         'ProdCOfromCH4                 ', 
                         'ProdCOfromNMVOC               ', 
                         'ProdCO2fromCO                 ',  
::
#==============================================================================
# %%%%% THE CloudConvFlux COLLECTION %%%%%
#
# Cloud convective flux (default = advected species)
#
# Available for all simulations
#==============================================================================
  CloudConvFlux.template:     '%y4%m2%d2_%h2%n2z.nc4',
  CloudConvFlux.frequency:    ${RUNDIR_HIST_TIME_AVG_FREQ}
  CloudConvFlux.duration:     ${RUNDIR_HIST_TIME_AVG_DUR}
  CloudConvFlux.mode:         'time-averaged'
  CloudConvFlux.fields:       'CloudConvFlux_?ADV?           ',
::
#==============================================================================
# %%%%% THE ConcAfterChem COLLECTION %%%%%
#
# Concentrations of OH, HO2, O1D, O3P immediately after exiting the KPP solver
# or OH after the CH4 specialty-simulation chemistry routine.
#
# OH:       Available for all full-chemistry simulations and CH4 specialty sim
#==============================================================================
  ConcAfterChem.template:     '%y4%m2%d2_%h2%n2z.nc4',
  ConcAfterChem.frequency:    ${RUNDIR_HIST_TIME_AVG_FREQ}
  ConcAfterChem.duration:     ${RUNDIR_HIST_TIME_AVG_DUR}
  ConcAfterChem.mode:         'time-averaged'
  ConcAfterChem.fields:       'OHconcAfterChem               ',
::
#==============================================================================
# %%%%% THE LevelEdgeDiags COLLECTION %%%%%
#
# Diagnostics that are defined on grid box level edges
#
# Available for all simulations
#==============================================================================
  LevelEdgeDiags.template:    '%y4%m2%d2_%h2%n2z.nc4',
  LevelEdgeDiags.frequency:   ${RUNDIR_HIST_TIME_AVG_FREQ}
  LevelEdgeDiags.duration:    ${RUNDIR_HIST_TIME_AVG_DUR}
  LevelEdgeDiags.mode:        'time-averaged'
  LevelEdgeDiags.fields:      'Met_CMFMC                     ',
                              'Met_PEDGE                     ',
                              'Met_PEDGEDRY                  ',
                              'Met_PFICU                     ',
                              'Met_PFILSAN                   ',
                              'Met_PFLCU                     ',
                              'Met_PFLLSAN                   ',
::
#==============================================================================
# %%%%% The StateMet COLLECTION %%%%%
#
# Fields of the State_Met object (also see the LevelEdgeDiags collection)
#
# Available for all simulations
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
                              'Met_FRLANDIC                  ',
                              'Met_FROCEAN                   ',
                              'Met_FRSEAICE                  ',
                              'Met_FRSNO                     ',
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
# %%%%% THE BoundaryConditions COLLECTION %%%%%
#
# GEOS-Chem boundary conditions for use in nested grid simulations
#
# Available for all simulations
#==============================================================================
  BoundaryConditions.template:   '%y4%m2%d2_%h2%n2z.nc4',
  BoundaryConditions.frequency:  00000000 030000
  BoundaryConditions.duration:   00000001 000000
  BoundaryConditions.mode:       'instantaneous'
  BoundaryConditions.fields:     'SpeciesBC_?ADV?             ',
::