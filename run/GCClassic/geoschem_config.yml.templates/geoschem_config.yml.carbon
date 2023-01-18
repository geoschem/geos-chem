---
### geoschem_config.yml: GEOS-Chem Runtime configuration options.
### Customized for simulations using the KPP carbon mechanism.

#============================================================================
# Simulation settings
#============================================================================
simulation:
  name: ${RUNDIR_SIM_NAME}
  start_date: [${RUNDIR_SIM_START_DATE}, ${RUNDIR_SIM_START_TIME}]
  end_date: [${RUNDIR_SIM_END_DATE}, ${RUNDIR_SIM_END_TIME}]
  root_data_dir: ${RUNDIR_DATA_ROOT}
  met_field: ${RUNDIR_MET}
  species_database_file: ./species_database.yml
  debug_printout: false
  use_gcclassic_timers: ${RUNDIR_USE_GCCLASSIC_TIMERS}

#============================================================================
# Grid settings
#============================================================================
grid:
  resolution: ${RUNDIR_GRID_RES_LONG}
  number_of_levels: ${RUNDIR_GRID_NLEV}
  longitude:
    range: ${RUNDIR_GRID_LON_RANGE}
    center_at_180: ${RUNDIR_CENTER_LON_180}
  latitude:
    range: ${RUNDIR_GRID_LAT_RANGE}
    half_size_polar_boxes: ${RUNDIR_GRID_HALF_POLAR}
  nested_grid_simulation:
    activate: ${RUNDIR_GRID_NESTED_SIM}
    buffer_zone_NSEW: ${RUNDIR_GRID_BUFFER_ZONE}

#============================================================================
# Timesteps settings
#============================================================================
timesteps:
  transport_timestep_in_s: ${RUNDIR_TRANSPORT_TS}
  chemistry_timestep_in_s: ${RUNDIR_CHEMISTRY_TS}

#============================================================================
# Settings for GEOS-Chem operations
#============================================================================
operations:

  chemistry:
    activate: true

  convection:
    activate: true

  pbl_mixing:
     activate: true
     use_non_local_pbl: ${RUNDIR_USE_NLPBL}

  transport:
    gcclassic_tpcore:                 # GEOS-Chem Classic only
      activate: true                  # GEOS-Chem Classic only
      fill_negative_values: true      # GEOS-Chem Classic only
      iord_jord_kord: [3, 3, 7]       # GEOS-Chem Classic only
    transported_species:
      - CH4
      - CO
      - CO2
      - OCS

#============================================================================
# Options for CH4
#============================================================================
CH4_simulation_options:

  use_observational_operators:
    AIRS: false
    GOSAT: false
    TCCON: false

  analytical_inversion:
    activate: false
    emission_perturbation: 1.0
    state_vector_element_number: 0
    use_emission_scale_factor: false
    use_OH_scale_factors: false

#============================================================================
# Options for CO
#============================================================================

CO_simulation_options:
  use_fullchem_PCO_from_CH4: true
  use_fullchem_PCO_from_NMVOC: true

#============================================================================
# Options for CO2
#============================================================================

CO2_simulation_options:

  sources:
    fossil_fuel_emissions: true
    ocean_exchange: true
    balanced_biosphere_exchange: true
    net_terrestrial_exchange: true
    ship_emissions: true
    aviation_emissions: true
    3D_chemical_oxidation_source: true

  tagged_species:
    save_fossil_fuel_in_background: false
    tag_bio_and_ocean_CO2: false
    tag_land_fossil_fuel_CO2: 
    tag_global_ship_CO2: false
    tag_global_aircraft_CO2: false
