# Changelog

This file documents all notable changes to the GEOS-Chem repository since version 13.4.1, including all GEOS-Chem Classic and GCHP run directory updates.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased 14.1.0]
### Added
  - Added dry deposition updates to Hg0 from Feinberg22 ESPI publication +
    AMAP emissions
  - Added MO2 + NO3 = NO2 + CH2O + HO2 reaction
  - Added capability to write species metadata to YAML file
  - Added satellite diagnostic (SatDiagn) collection, to archive several
    fields within a user-defined local-time interval. CAVEAT: For now,
    only one local-time interval is permitted.
  - Added GCHP run script and environment files for MIT clusters Hex and Svante

### Changed
  - Moved in-module variables in global_ch4_mod.F90 to State_Chm
  - Moved in-module variables in hco_interface_gc_mod.F90 to State_Met and State_Chm

### Fixed
  - Fixed sign of Arrhenius "A" coefficient in reaction ETO = HO2 + 2CH2O
  - Fixed products in HOBr + SO2 and HOCl + SO2 reactions
  - Changed MW_g value of CH4 from 16.05 to 16.04
  - Added "WD_CoarseAer:true" for SO4s and NITs in species_database.yml
  -Fixed bug in computing State_Met surface type logicals (IsLand, IsWater, etc)


## [14.0.2] - 2022-11-29
### Fixed
  - Added fix for writing dry-run header to log file
  - Updated KPP diagnostics archive flags
  - Rewrote code to avoid memory leaks (identified by the code sanitizer)
  - Updated EDGAR v6 CH4 emission files to correct timestamp issue
  - Updated CH4 Lakes emission files to correct time unit issue
  - Added fix for CH4_RICE emissions from EDGAR v6
  - Fixed indentation error in the `legacy_bpch` section of
    `geoschem_config.yml` template files
  - Removed "dry air" from the metadata of fields `State_Met%AIRVOL` and
    `State_Met%BXHEIGHT`
  - Applied fixes for CESM runs: Turned off sea salt emissions; Modified time
    cycle flag for YUAN_MODIS_LAI

### Changed
  - Updated CESM HISTORY.rc to work with new CESM-GC diagnostics interface
  - Updated sample fullchem restart files copied to run directories to 14.0.0
    10-year benchmark output


## [14.0.1] - 2022-10-31
### Fixed
  - Corrected units in metadata for State_Met%AirNumDen and State_Met%PHIS
  - Fixed file path for AEIC2019_DAILY emissions for aerosol-only simulations
  - Fixed GCHP bug to populate non-species data in mid-run restart files
  - Fixed typo preventing ND51 satellite diagnostic from turning on

### Changed
  - Documented and cleaned up GCHP run script operational examples
  - Updated README.md and AUTHORS.txt
  - Set species concentration arrays as pointers to internal state in GCHP
  - Updated Restart collection in HISTORY.rc to save out BXHEIGHT and TROPLEV for all simulations


## [14.0.0] - 2022-10-25
### Added
  - Added user registration with dynamodb database during run directory creation
  - Added Hg simulation with KPP
  - Added yaml-format config file geoschem_config.yml which replaces input.geos
  - Added native GEOS-FP and mass fluxes options to GCHP run directory creation
  - Added cap_restart file to GCHP run directories to set simulation start time
  - Added updates for compatibility with CESM, GEOS, and WRF-GC

### Fixed
  - Fixed missing output boundary conditions on first timestep of run
  - Added missing entries for POG1, POG2, and pFe to HEMCO_Config.rc
  - Reverted GC-Classic pressure fixer to v13.3 to fix bug in v13.4
  - Fixed dry deposition of methanol over oceans
  - Fixed issues in creating run directory for GCAP2
  - Removed duplicate species for SO4 in aciduptake.eqn
  - Fixed CEDS_CO2_SHP emissions in HEMCO_Config.rc file for CO2 simulation
  - Fixed Volcano_Table entry in HEMCO config template for GCHP
  - Fixed transport tracers simulation in GCHP
  - Applied fix to avoid divide-by-zero in routine MMR_Compute_FLux
  - Fixed HEMCO diagnostic counter zero warnings in full chemistry simulation
  - Fixed bug in totalOC diagnostic
  - Fixed bugs causing differences when splitting up GC-Classic and GCHP simulations in time
  - Fixed bug setting GEOS-FP meteorology in GCHP run directories

### Changed
  - Updated KPP to version 2.5.0
  - Updated GCHP run scripts to easily segment runs in time
  - Changed GCHP restart filename convention to exclude seconds
  - Updated offline biogenic VOC and soil NOx emissions
  - Reduced root logging level for MAPL from INFO to WARNING
  - Changed 4D State_Chm%Species array to vector of 3D concentration arrays
  - Renamed GCHP config file runConfig.sh to setCommonRunSettings.sh
  - Moved restart file location in run directory to Restarts subdirectory
  - Updated sample restart files copied to run directories to 14.0.0
    1-year benchmark output

### Removed
  - Removed TMPU1, SPHU1, PS1_WET, and PS1_DRY from GC-Classic restart file
  - Removed input.geos; replaced with geoschem_config.yml
  - Removed HEMCO.log output file; HEMCO log info now sent to main GEOS-Chem log
