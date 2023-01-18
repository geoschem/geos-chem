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
  debug_printout: false

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

  legacy_bpch:                #          1         2         3
     output_menu:             # 1234567890123456789012345678901
       schedule_output_for_JAN: 3000000000000000000000000000000
       schedule_output_for_FEB: 30000000000000000000000000000
       schedule_output_for_MAR: 3000000000000000000000000000000
       schedule_output_for_APR: 300000000000000000000000000000
       schedule_output_for_MAY: 3000000000000000000000000000000
       schedule_output_for_JUN: 300000000000000000000000000000
       schedule_output_for_JUL: 3000000000000000000000000000000
       schedule_output_for_AUG: 3000000000000000000000000000000
       schedule_output_for_SEP: 300000000000000000000000000000
       schedule_output_for_OCT: 3000000000000000000000000000000
       schedule_output_for_NOV: 300000000000000000000000000000
       schedule_output_for_DEC: 3000000000000000000000000000000

    gamap:
      diaginfo_dat_file: ./diaginfo.dat
      tracerinfo_dat_file: ./tracerinfo.dat

    bpch_diagnostics:
      TOMAS_aerosol_emissions: "0 all"
      TOMAS_rate: "0 all"
      TOMAS_3D_rate: "0 all"
      ND65_prodloss:
        activate: true
        number_of_levels: 72

    ND51_satellite:
      activate: false
      output_file: ts_satellite.YYYYMMDD.bpch
      tracers:
        - 1
        - 2
        - 501
      UTC_hour_for_write: 0
      averaging_period_in_LT: [9, 11]
      IMIN_and_IMAX_of_region: [1, 72]
      JMIN_and_JMAX_of_region: [1, 46]
      LMIN_and_LMAX_of_region: [1, 1]

    ND51b_satellite:
      activate: false
      output_file: ts_13_15_NA..YYYYMMDD.bpch
      tracers:
        - 1
        - 2
        - 501
      UTC_hour_for_write: 1
      averaging_period_in_LT: [13, 15]
      IMIN_and_IMAX_of_region: [1, 72]
      JMIN_and_JMAX_of_region: [1, 46]
      LMIN_and_LMAX_of_region: [1, 1]
