# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased 14.0.0]
### Added
- Register first time users with dynamodb database feature
- Add Hg simulation with KPP
- Save boundary conditions on first timestep
- Convert input.geos to YAML format; new file is geoschem_config.yml
- Add native GEOS-FP and mass fluxes options in GCHP run directory creation
- Add GEOS-Chem updates for GEOS from GMAO
- Add GEOS-Chem updates for CESM

### Fixed
- Specify AOD wavelength 999 nm for use in Fast-JX
- Add missing entries for POG1, POG2, and pFe to HEMCO_Config.rc
- Revert GC-Classic pressure fixer to v13.3
- Add special treatment for MOH in dry deposition if not over land
- Fix issues in creating run directory for GCAP2
- Remove duplicate species for SO4 in aciduptake.eqn bug
- Update CEDS_CO2_SHP emissions in HEMCO_Config.rc file for CO2 simulation
- Fix Volcano_Table in HEMCO config template for GCHP
- Fix transport tracers simulation in GCHP
- Avoid divide-by-zerio in routine MMR_Compute_FLux
- Fix HEMCO diagnostic counter zero warnings in full chemistry simulation
- Fix bug in totalOC diagnostic
- Fix bugs to remove GCHP diffs if splitting up simulations in time

### Changed
- Update GCHP run directory files to easily segment runs in time
- Update to KPP 2.5.0
- Change GCHP restart filename convention to exclude seconds
- Update offline biogenic VOC and soil NOx emissions
- Reduce root logging level for MAPL to WARNING
- Change 4D State_Chm%Species array to vector of 3D concentration arrays feature

### Removed
- Remove TMPU1, SPHU1, PS1_WET, and PS1_DRY from GC-Classic restart file
