---
### geoschem_config.yml: GEOS-Chem Runtime configuration options.
### Customized for simulations using the KPP carbon mechanism.

#============================================================================
# Simulation settings
#============================================================================
simulation:
  name: ${RUNDIR_SIM_NAME}
  chem_inputs_dir: ./ChemDir/
  species_database_file: ./species_database.yml
  species_metadata_output_file: OutputDir/geoschem_species_metadata.yml
  verbose:
    activate: false
    on_cores: root       # Allowed values: root all

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
    state_vector_element_number: 0
    emission_perturbation_factor: 1.0
    OH_perturbation_factor: 1.0
    perturb_CH4_boundary_conditions: false
    CH4_boundary_condition_ppb_increase_NSEW: [0.0, 0.0, 0.0, 0.0]
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

#============================================================================
# Settings for diagnostics (other than HISTORY and HEMCO)
#============================================================================
extra_diagnostics:

  obspack:
    activate: false
    quiet_logfile_output: false
    input_file: ./obspack_co2_1_OCO2MIP_2018-11-28.YYYYMMDD.nc
    output_file: ./OutputDir/GEOSChem.ObsPack.YYYYMMDD_hhmmz.nc4
    output_species:
      - CH4
      - CO
      - CO2

  planeflight:
    activate: false
    flight_track_file: Planeflight.dat.YYYYMMDD
    output_file: plane.log.YYYYMMDD