# Changelog

This file documents all notable changes to the GEOS-Chem repository starting in version 14.0.0, including all GEOS-Chem Classic and GCHP run directory updates.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [14.5.2] - 2025-02-12
### Added
- Implemented the Global Rice Patty Inventory (GRPI) for CH4 and carbon simulations to replace EDGAR rice emissions
- Added run directory creation for processed cubed-sphere GEOS-IT meteorology
- Added GC-Classic and GCHP environment files, build scripts, and run scripts for MSU Orion cluster

### Changed
- Updated GC-Classic and GCHP environment files, build scripts, and run scripts for NASA discover cluster
- Updated rundir scripts to ask for confirmation before building the KPP-Standalone executable
- Updated rundir scripts to print a reminder to compile with `-DKPPSA=y` to build the KPP-Standalone executable
- Updated `integrationTestCreate.sh` and `parallelTestCreate.sh` scripts to decline building the KPP-Standalone.

### Fixed
- Fixed GCHP refresh time for `CO2_WEEKLY` scale factors so updated daily
- Fixed bug in GCHP GEOS-IT run directory using raw lat-lon fields on NASA discover cluster

## [14.5.1] - 2025-01-10
### Added
- Added Australian Hg emissions for 2000-2019 from MacFarlane et. al. [2022], plus corresponding mask file
- Added comments in GEOS-Chem Classic `HISTORY.rc` template files advising users not to change the `BoundaryConditions.frequency` setting
- Added `.zenodo.json` for auto-DOI generation upon version releases

### Fixed
- Reverted CH4 livestock emissions to EDGAR v7 to avoid hotspots and to apply seasonality

### Removed
- Removed unused RUNDIR settings for GCHP pressure units and scaling

## [14.5.1] - 2025-01-10
### Added
- Added allocate guards for arrays in `pressure_mod`
- Added `State_Diag%SatDiagnEdgeCount` counter for the `SatDiagnEdge` collection
- Added `State_Diag%Archive_SatDiagnEdgeCount` field
- Added `State_Diag%Archive_SatDiagnEdge` field
- Added routine `SatDiagn_or_SatDiagnEdge` in `History/history_utils_mod.F90`
- Added error trap in `History/history_mod.F90` to ensure that collection duration is always shorter than frequency
- Added KPP standalone interface (archives model state to selected locations)
- Added `https://github/geoschem/KPP-Standalone` as a Git submodule
- Added comments in `./run/sharedcleanRunDir.sh` describing the `--force` option (i.e. remove files w/o user confirmation)
- Specified meteorology source in GCHP geoschem_config.yml
- Added Input_Opt logical for whether to reconstruct convective precipitation fluxes rather than use met-fields
- Added to run directory creation a warning about convection discontinuity and bug if GEOS-FP meteorology is chosen
- Added surface precipitation flux fields as inputs to GCHP

### Changed
- Renamed `Emiss_Carbon_Gases` to `CO2_Production` in `carbon_gases_mod.F90`
- Updated start date and restart file for CO2 and tagCO simulations for consistency with carbon simulations
- Allocated `State_Diag%SatDiagnPEDGE` ffield with vertical dimension `State_Grid%NZ+1`
- Modified `run/GCClassic/cleanRunDir.sh` to skip removing bpch files, as well as now removing `fort.*` and `OutputDir/*.txt` files
- Edited `run/shared/kpp_standalone_interface.yml` to include additional entries under `active cells` and `locations`
- Changed doing Linoz and Linearized chemistry messages to print only if verbose
- Updated HEMCO subroutine calls for error and log handling changes in HEMCO 3.9.1
- Updated configuration files for using GEOS-Chem 14.5 in CESM
- Modified tagCO simulation to use GFED4 biomass burning emissions and GEOS-Chem v5 OH fields for consistency with carbon simulation
- Changed integration tests to use Harvard Cannon GNU 12 environment files by default

### Fixed
- Added a fix to skip the call to KPP when only CO2 is defined in the carbon simulation
- Added fix to turn on ship emissions for CO2 in the carbon simulation
- Updated `HEMCO_Config.rc` for carbon simulation to read data based on carbon species used
- Fixed entries for CO2 emissions in `ExtData.rc.carbon`
- Fixed metals simulation name in config file template comments
- Fixed bug in `download_data.py` which caused script to fail if log filename contained uppercase characters.
- Fixed the satellite diagnostics counters from being inadvertently being reset
- Fixed segmentation fault in qfyaml when running with certain compilers without debug flags on
- Fixed errors in adjoint-only code preventing successful adjoint build
- Fixed zero convective precipitation and high cloud base in runs using GEOS-FP (>=01Jun2020) or GEOS-IT
- Updated GEOS-only code and configuration files for compatibility with GEOS-Chem 14.5
- Fixed missing Is_Advected for TMB in species_database.yml
- Fixed typos in `HEMCO_Config.rc` for CH4 simulations causing mobile combustion emissions to be double counted
- Fixed handling of FIRST flag in carbon_gases_mod.F to limit log prints to first timestep only
- Removed extraneous pressure correction in GCHP carbon simulations by adding 'activate: true' to geoschem_config.yml
- Fixed bug in GC-Classic OCS emissions where unit conversion of km2 to m2 occurred twice
- Changed dimension of EmisOCS_Total from 3D to 2D since all emissions for all sectors are 2D
- Added fixes to only apply archived PCO_CH4 field for carbon simulations with CO only

### Removed
- Removed duplicate `WD_RetFactor` tag for HgClHO2 in `species_database.yml`
- Removed error messages in HEMCO interface pointing users to HEMCO log

## [14.5.0] - 2024-11-07
### Added
- Added vectors `State_Chm%KPP_AbsTol` and `State_Chm%KPP_RelTol`
- Added setting `KPP_AbsTol` to 1e5 for dummy species in `species_database.yml` and `species_database_hg.yml`
- Implemented PPN photolysis from Horner et al (2024)
- Added four new species ALK4N1, ALK4N2, ALK4O2, and ALK4P to address issues in ALK4 and R4N2 chemistry following Brewer et al. (2023, JGR)
- Added new species ALK4N1 and ALK4N2 to Ox family in KPP
- Added Cloud-J input parameters to geoschem_config.yml in new photolysis sub-menu called cloud-j
- Added computation of water concentration to use in photolysis for application of UV absorption by water in Cloud-J v8
- Added ACO3, ACR, ACRO2, ALK4N{1,2,O}2, ALK4P, ALK6, APAN, APINN, APINO2, APINP, AROCMCHO, AROMCO3, AROMPN, BPINN, BPINO2, BPINON, BPINOO2, BPINOOH, BPINP, BUTN, BUTO2, C4H6, C96N, C96O2, C9602H, EBZ, GCO3, HACTA, LIMAL, LIMKB, LIMKET, LIMKO2, LIMN, LIMNB, LIMO2H, LIMO3, LIMO3H, LIMPAN, MEKCO3, MEKPN, MYRCO, PHAN, PIN, PINAL, PINO3, PINONIC, PINPAN, R7N{1,2}, R7O2, R7P, RNO3, STYR, TLFUO2, TLFUONE, TMB, ZRO2 to `species_database.yml` following Travis et  al. 2024
- Added TSOIL1 field to `State_Met` for use in HEMCO soil NOx extension. This should only be read in when the `UseSoilTemperature` option is true in HEMCO config

### Changed
- Copied values from `State_Chm%KPP_AbsTol` to `ATOL` and `State_Chm%KPP_RelTol` to `RTOL` for fullchem and Hg simulations
- Introduced seasalt Ca, K, Mg back to aerosol thermodynamics via HETP.
- Updated `HEMCO_Config.rc.fullchem` (GCClassic + GCHP) and `ExtData.rc` to add emissons of new species from Travis et al 2023
- Activated the `DryDep` collection for GCClassic & GCHP fullchem benchmarks
- Reduced the GCHP `DryDep` collection to only the necessary species for benchmarks
- Removed unused `VDIFFAR` routine from `vdiff_mod.F90`
- Updated MW for CH4 and OH in `global_ch4_mod.F90`
- Added fix to not convert from kg/kg to mol/mol before passing State_Chm to PBL mixing in `vdiff_mod.F90`
- Updated GC-Classic and GCHP run scripts and environment files for NASA discover cluster
- Updated `GFED4_Climatology` entries to point to the climatology file for 2010-2023
- Moved aerosol optical properties files to a new data directory specified in geoschem_config.yml rather than specifying in photolysis input files
- Moved calls to `RD_AOD` and `CALC_AOD` from `Init_Aerosol` rather than `Init_Photolysis`
- Updated ResME CH4 reservoir emissions to apply seasonality via mask file
- Changed fullchem restart file folder from `GC_14.3.0` to `GC_14.5.0`
- Excluded HEMCO interface and ExtState fields from `MODEL_CESM` in `hco_interface_gc_mod.F90` for compatibility with CESM, which runs HEMCO separately

### Fixed
- Simplified SOA representations and fixed related AOD and TotalOA/OC calculations in benchmark
- Changed mass conservation adjustment in `vdiff_mod.F90` to use a mass tendency with units of `kg species/kg dry air`
- Converted the top pressure edge from hPa to Pa in `vdiff_mod.F90`
- Updated `Jval_` entries in `run/GCHP/HISTORY.rc.templates/HISTORY.rc.fullchem`
- Updated species database Is_Photolysis entries to remove J-value diagnostics with all zeros in full chemistry simulation
- Removed EDGAR8_CH4_AWB emissions from CH4 and carbon simulations to avoid double counting with GFED
- Fixed formatting error in `.github/workflows/stale.yml` that caused the Mark Stale Issues action not to run
- Fixed emissions in GCHP carbon ExtData.rc so that data in molecules/cm2/s are converted to kg/m2/s

### Removed
- Removed dry-run checks for files that are no longer needed for Cloud-J v8 from `cldj_interface_mod.F90`

## [14.4.3] - 2024-08-13
### Added
- Added tropopause pressure field in the satellite diagnostic (by @eamarais)
- Added ModelEe.2 (GCAP 2.0) simulation to integration tests for GCClassic
- Added simulation with all diagnostics on in HISTORY.rc to integration tests for GCClassic (including Planeflight + ObsPack) and GCHP
- Added descriptive error message in `Interfaces/GCHP/gchp_historyexportsmod.F90`
- Auto-update GCHP HEMCO_Diagn.rc settings at run-time to ensure seasalt, dust, soil NOx, and biogenic emissions match settings in HEMCO_Config.rc

### Fixed
- Added brackets around `exempt-issue-labels` list in `.github/workflows/stale.yml`

### Removed
- Removed `XNUMOL_H2O2 / CM3PERM3` in routine `Chem_H2O2`, which removes an unnecessary unit conversion for the aerosol-only simulation

## [14.4.2] - 2024-07-24
### Added
- Added number of levels with clouds for photolysis to geoschem_config.yml and Input_Opt to pass to Cloud-J
- Added `State_Grid%CPU_Subdomain_ID` and `State_Grid%CPU_Subdomain_FirstID` as "identifier numbers" for multiple instances of GEOS-Chem on one core in WRF and CESM
- Added transport tracer run directory option for global half-degree GC-Classic run with GEOS-IT 0.5x0.625 fields

### Changed
- Now reset `State_Diag%SatDiagnCount` to zero in routine`History_Write` (instead of in `History_Netcdf_Write`)
- Update rundir creation scripts to turn off the MEGAN extension for "standard" fullchem simulations
- Updated emissions used in CESM to match standard emissions used in the 14.4 offline model
- Disable support For FAST-JX for all simulations except Hg
- Only read photolysis data in `Init_Photolysis` in first instance of GEOS-Chem on each PET in CESM as PIO requires it
- Replace calls to `GEOS_CHEM_STOP` with calls to `GC_Error` in `planeflight_mod.F90`
- Script `test/integration/GCHP/integrationTestExecute.sh` now resets `cap_restart` time to `000000`, to facilitate manual restart

### Fixed
- In `Headers/roundoff_mod.F90`, first cast and then only round off if `places > 0`
- Typo in `setCommonRunSettings.sh` that made GCHP always choose mass fluxes for meteorology
- Fixed bug in # levels with cloud used in photolysis when using GCAP met or CESM
- Fixed typos for `SatDiagnEdge` collection in `HISTORY.rc` templates
- The `SatDiagnOH` diagnostic now works for the carbon simulation
- Restored missing fields for `UVFlux` collection in `run/GCClassic/HISTORY.rc.templates/HISTORY.rc.fullchem`
- Comment out `UVFlux` diagnostic in the "alldiags" integration test, there is a floating point error.  Look at this later.
- Now use SO4 instead of O3 in the GCHP fullchem budget diagnostic (SO4 is soluble, O3 is not)
- Convert `UVFlux_Tag_Names` to uppercase in the comparison in `Get_UVFlux_Bin`(located in`Headers/state_diag_mod.F90`)
- Fixed typo (missing `_` character) in GCHP `DryDep` collection diagnostic entries
- Commented out with `###` emissions diagnostics in the GCHP `HISTORY.rc.fullchem` template that are not present in the corresponding `HEMCO_DIAGN.rc` template

### Removed
- Entry `SatDiagnPEDGE` from the `SatDiagn` collection; This needs to go into the `SatDiagnEdge` collection.

## [14.4.1] - 2024-06-28
### Added
- Added initialization of PHOTDELTA in `ucx_h2so4phot` to avoid run-time error in CESM
- Added Cloud-J status output and error handling for it

### Changed
- Alphabetically sort Complex SOA species into `geoschem_config.yml` in run directory creation 
- Use hard-coded years for met fields and BC files in `HEMCO_Config.rc` so they are not read hourly
- Updated `run/CESM` with alphabetical sorting of species in `geoschem_config.yml`
- Added clarifying comments in GCHP configuration files for several settings, particularly related to domain decomposition, mass fluxes, and stretched grid
- Added pre-run GCHP configuration checks to `setCommonRunSettings.sh` related to domain decomposition, mass fluxes, and stretched grid.
- Changed search criteria for GCHP auto-update of met-field refresh frequency to not rely on presence of `MetDir` symlink in `ExtData.rc` file path

### Fixed
- Fixed formatting error in `.github/workflows/stale.yml` that caused the Mark Stale Issues action not to run
- Fixed typo `$GCAPVERTRESL` -> `$GCAPVERTRES` in `HEMCO_Config.rc.fullchem` template file
- Fixed GCHP `ExtData.rc` entry for lightning climatology files

### Removed
- Removed `BudgetWetDep*` entries from simulations with no soluble species in `HISTORY.rc` templates
- Disabled `run/CESM` ParaNOx extension by default in `HEMCO_Config.rc`
- Removed MPI broadcasts in CESM-only UCX code; MPI broadcast done at coupler level
- Remove enabling O-server in GCHP for high core counts

### Fixed
- In `Headers/roundoff_mod.F90`, first cast and then only round off if `places > 0`

## [14.4.0] - 2024-05-30
### Added
- Added `SpcConc%Units` for species-specific unit conversion
- Diel and day-of-week scale factors for CEDS global base emissions
- `Input_Opt%Satellite_CH4_Columns` logical flag; Set this to true if any of AIRS, GOSAT, TCCON observational operators are selected
- Add explicit handling of gravitational settling and hygroscopic growth in dry deposition
- Added CO2, CO, and OCS single-tracer carbon simulations to the integration tests
- Added missing entry in `HEMCO_Config.rc` for natural gas postmeter CH4 emissions in GHGIv2 Express Extension
- Added tagged species capability and PM25nit and PM25nh4 diagnostics for GEOS runs
- Added `real*4` diagnostics for State_Met logical masks IsWater, IsLand, IsIce, and IsSnow
- New parameterization for effective radius of SNA/OM aersols (see PR #2236)
- New `CHEM_INPUTS/FAST_JX/v2024-05` and `CHEM_INPUTS/FAST_JX/v2024-05-Hg` folders with updated `org.dat` and `so4.dat` files
- Added global continental chlorine (pCl and HCl) emissions
- Extended GFED4 emissions through the end of 2023
- Added a parameterization for dry aerosol size (Rg) for SNA and OM aerosols. Updated AOD calculation reflecting varying aerosol size.

### Changed
- Updated routines in `GeosUtil/unitconv_mod.F90` for species-specific unit conversion
- Halt timers during calls to `Convert_Spc_Units` so as to time unit conversions separately
- Streamline `IF` statements for CH4 observational operators in `Interfaces/GCClassic/main.F90`
- Disable parallel loop in `Do_Convection` when using TOMAS; it causes unit conversion issues.  Revisit later.
- Add explicit handling of gravitational settling and hygroscopic growth in dry deposition
- Added CO2, CO, and OCS single-tracer carbon simulations to the integration tests
- GitHub Action config file `.github/workflows/stale.yml`, which replaces StaleBot
- Switch from fixed to monthly timezones, which account for daylight savings time more accurately when computing emissions
- Updated NOAA GMD surface CH4 boundary conditions through 2022
- Rename `NITs_Jscale_JHNO3` to `NITs_Jscale` and `NIT_Jscale_JHNO2` to `NIT_Jscale` in `geoschem_config.yml` templates
- Updated volcano emissions from GMAO v202005 product to v202401 which extends to the end of 2024
- Use local scale height and level thickness to determine the PBL to determine the PBL top level and PBL pressure thickness
- Update drydep mean diameters of aerosols to account for size distribution
- Corrected the formula for 1st order heterogeneous chemical loss on stratospheric aerosol for NO2, NO3, and VOC.
- Fixed incorrect time refresh entries and other errors in `run/GCHP/ExtData.rc.templates/ExtData.rc.carbon`
- Changed time range entries in HEMCO_Config.rc for met, restart, and BC files to use year, month, and day tokens instead of hardcoded range
- Renamed `State_Met%FRSNO` and `State_Met%FRLANDIC` to `State_Met%FRSNOW` and `State_Met%FRLANDICE`
- Renamed isorropiaII_mod.F90 to aerosol_thermodynamics_mod.F90
- Changed aerosol thermodynamics scheme from ISORROPIA II to HETP for fullchem and APM
- Changed input data paths in `run/GEOS` directory to match location change on NASA discover cluster
- Use new mask files at 0.1 x 0.1 degree resoluiton for CH4/tagCH4/carbon simulations to avoid I/O bottlenecks
- Update config files for CH4/carbon simulations to avoid reading the same variable multiple times
- Converted Github issue templates to issue forms using YAML definition files

### Fixed
- Corrected the formula for 1st order heterogeneous chemical loss on stratospheric aerosol for NO2, NO3, and VOC.
- Use rate-law function `GCARR_ac` for rxns that have Arrhenius `B` parameters that are zero
- Now use correct index `WEAEROSOL(I,J,L,2+NDUST)` in routine `Settle_Strat_Aer` of `GeosCore/ucx_mod.F90`
- Now get density of BCPI species from the species database in `ucx_mod.F90`
- Fix issues that prevented single-species carbon simulations from running in GCHP
- Update `HEMCO_Config.rc.carbon` and `ExtData.rc.carbon` templates for consistency
- Updated several emissions files for CO and CH4 for COARDS and MAPL compliance
- Fixed several issues in GCHP single-species carbon simulation setup scripts
- Corrected the formula for 1st order heterogeneous chemical loss on stratospheric aerosol for NO2, NO3, and VOC.
- Corrected the formula for 1st order heterogeneous chemical loss on stratospheric aerosol for NO2, NO3, and VOC.
- Change restart file time cycle flag from `EFYO` to `CYS` for TOMAS simulations to avoid missing species error.
- Now define `REEVAPSO2` in wetscav_mod when units are kg species; this avoids floating-point errors.
- Fixed `State_Met%FRSNO` to be fraction of grid box with snow rather than fraction of land with snow
- Fixed variable definitions in the `DryDep` collection of `run/GCHP/HISTORY.rc.templates/HISTORY.rc.fullchem`

### Removed
- Legacy binary punch diagnostic code contained within `#ifdef BPCH_DIAG` blocks
- `IU_BPCH` logical file unit (in `GeosUtil/file_mod.F90`)
- Removed tagged CH4 and CO species handling from `carbon_gases_mod.F90`
- GitHub config files `.github/stale.yml` and `.github/no-response.yml`
- Unused CO2 and carbon simulation options from `geoschem_config.yml` (and from related code in co2_mod.F90).
- Removed ISORROPIA
- Removed `Begin` array in do_fullchem (declared but not used)
- Removed tagCH4 simulation as option
- Removed `--request-payer requester` from `run/shared/download_data.py`; the `s3://gcgrid` data is open-source

## [14.3.1] - 2024-04-02
### Added
- Added operational run scripts for the Imperial College London (ICL) cluster
- Added new vertical region option to budget diagnostic for fixed bottom and top levels
- Added GEOS-IT processed lat-lon fields as a valid option when creating GCHP run directories
- Functions `charArr2str` and `str2CharArr` in `Headers/charpak_mod.F90`
- Field `State_Diag%Obspack_CharArray` as a 2-D character array
- Added util folder in run/CESM to include .cdl file used to generate CESM NetCDF input file for deposition
- Add GCClassic operational example environment files for Harvard Cannon
- Added new GCHP history collections for advection diagnostics
- Added slash in front of names of LUT files read into `photolysis_mod.F90` to avoid needing it in path

### Changed
- Updated Harvard Cannon operational run scripts to use `huce_cascade` instead of `huce_intel`; also added `sapphire`
- Changed exponent 'e' to 'd' for one entry in KPP to prevent precision error in external models
- Changed GCHP sample run scripts to not print script execution commands to log
- Changed offline emissions grid resolution templates in config files to be more descriptive
- Read `obspack_id` from netCDF files into a character array, then convert to string
- Add `#SBATCH -c 1` to GCHP integration test scripts and sample run scripts for Harvard Cannon
- In GCC/GCHP integration tests, passing `-s none` will run compile-only tests.  Query user to proceed or to exit.
- GCC/GCHP integration tests will exit immediately if `scheduler` is omitted.
- Now use `raw` instead of `native` in GCHP run directory scripts & templates
- Rename env var `RUNDIR_METLIGHTNING_DIR_NATIVE` to `RUNDIR_METLIGHTNING_DIR`
- Rename env var `RUNDIR_METLIGHTNING_NATIVE_RES` to `RUNDIR_METLIGHTNING_RES`
- Updated config files used in CESM from GEOS-Chem 14.1 to 14.3
- Don't create run directories for integration/parallel tests if invoked with `-t compile`
- Refactor integration and parallel test scripts to reduce the number of input arguments
- Copy utility scripts that allow you to resubmit failed to integration and parallel test root directories
- Update GCHP operational example environment files for Harvard Cannon
- Do not run GCClassic integration test compile jobs in the background
- Updated integration tests to pass quick option to compile scripts
- Removed emissions handling from `global_ch4_mod.F90` and `carbon_gases_mod.F90` and instead apply scale factors to emissions directly in `HEMCO_Config.rc`
- Loop over advected species CH4 chemistry routines to allow for multiple CH4 tracers within analytical inversion framework
- Updated CH4 global anthropogenic emission inventory from EDGARv7 to EDGARv8

### Fixed
- Fixed unit conversions in GEOS-only code
- Fixed GEOS-IT native lat-lon filenames used for clusters other than discover
- Fixed offline emission paths set when using GEOS-IT meteorology
- Fixed format issue in input_mod RRTMG print statement caught by some compilers
- Fixed GEOS-IT SLP and TROPP scaling in pre-processed files used in GCHP
- Fixed reading of NEI emissions through HEMCO
- Fixed incorrect units metadata for `State_Met%PHIS`
- Fixed bug in transport tracer ST80 mask criteria which prevented mask from ever being zero

### Removed
- Removed MPI broadcasts in CESM-only photolysis code; will read on all cores
- Removed `State_Chm%CH4_EMIS`

## [14.3.0] - 2024-02-07
### Added
- Added capability for TOMAS simulations in GCHP
- Added State_Chm%nTomasBins to replace hardcoded bins in TOMAS diagnostics
- Added interface to Cloud-J package for computing photolysis rates
- Added compile-time option FASTJX to use legacy Fast-JX photolysis instead of Cloud-J
- Added new diagnostics OD600 and TCOD600 for 600 nm optical depths (per-level and total column) used for computing J-values in either Fast-JX or Cloud-J
- Added GEOS-IT as meteorology option and labeled as beta during run directory creation until full inventory and offline emissions are available.
- Added support for running GEOS-Chem on the NASA discover cluster
- Added inclusion of c30 restart file in GCHP run directories since c24 and c48 not supported when using GEOS-IT meteorology
- Added automatic updating of GCHP lightning climatology in ExtData.rc based on settings in HEMCO_Config.rc
- Added two new diagnostics to track number of negative concentrations after first and last KPP integration
- Added capability of running GEOS-Chem transport tracer simulation within the GEOS model
- Added radiative forcing contributions due to trop-only ozone, CFCs, water vapor, N2O, CO2 and changes in stratosphere to RRTMG
- Added computation of radiative forcing at the tropopause to RRTMG
- Added option to compute stratospherically-adjusted radiative forcing at the tropopause using RK4 time marching integration with fixed dynamical heating approximation (FDH)
- Added experimental option to apply seasonally-evolving fixed dyanmical heating approximation in RRTMG

### Changed
- Updated fullchem mechanism following JPL/IUPAC. See `KPP/fullchem/CHANGELOG_fullchem.md` for details.
- Reorganized GCHP run directory creation prompts for GEOS-FP native meteorology input
- Converted TOMAS bpch diagnostics to netCDF
- Now read the Hg restart file from `ExtData/GEOSCHEM_RESTARTS/v2023-12`
- Increse requested time limits in GCHP integration tests (compile 2h30m, run 5h)
- Changed CO2 concentration used in RRTMG to be modifiable in geoschem_config.yml
- Changed water vapor used in RRTMG to match to tracer field at all altitudes
- Updated restart file path for GCHP TOMAS simulations
- Look for fullchem restarts in the `GEOSCHEM_RESTARTS/GC_14.3.0` folder
- Look for fullchem/aerosol boundary conditions in the `HEMCO/SAMPLE_BCs/GC_14.3.0/fullchem` folder

### Fixed
- Fixed bug in stratospheric aerosols optical depths passed to Fast-JX
- Restored consideration of both isSnow and isIce in dry deposition
- Fixed calculation of `FRLAND_NOSNO_NOICE` in `calc_met_mod.F90`
- Added missing units in comments of `KPP/fullchem/commonIncludeVars.H`
- Use run directory (not absolute path) to determine the executable file name in integration & parallel tests.
- Fixed memory leaks in `State_Chm%AerMass` and `State_Chm%Phot` containers
- Fixed incorrect time-avaging in RRTMG diagnostics wheres zeros included prior to first RRTMG call
- Added fix for runaway HMS chemistry. See `KPP/fullchem/CHANGELOG_fullchem.md` for details.

### Removed
- Removed references to unused met-fields RADLWG and LWGNT
- Removed inclusion of c360 restart file in GCHP run directories
- Reduced timers saved out to essential list used for benchmarking model performance
- Removed `State_Chm%Spc_Units`; this is now superseded by `State_Chm%Species(:)%Units`

## [14.2.3] - 2023-12-01
### Added
- GEOS-Chem Classic rundir script `run/GCClassic/setupForRestarts.sh`

### Changed
- Added the `-n` aka `--no-bootstrap` option to integration tests to disable bootstrapping missing species in restart files
- Use integer parameters for species units instead of strings (for computational efficiency)
- Update error message for missing surface CH4 emissions with instructions on how to resolve the problem
- Change GCHP grid resolution threshold for lowering timesteps from C180 inclusive to C180 exclusive
- Read GEOS-Chem Classic restart file paths from the relevant `download_data.yml` file
- Moved aerosol_mod module variables to new State_Chm container called AerMass

### Fixed
- Prevent `POAEMISS` from being assigned a value if not allocated (in `carbon_mod.F90`)
- Changed incorrect comment about static H2O option in `GeosCore/input_mod.F90`
- Fixed typos (`GCClassic` -> `GCHP`) written to GCHP integration test log files
- Add fix to properly read GHGI v2 express extension emissions in CH4 and carbon simulations
- Move OH perturbation scale factor to outside EMISSIONS logical bracket in HEMCO_Config.rc files for CH4 and carbon simulations

### Removed
- Remove definition of METDIR from primary `HEMCO_Config.rc` files to ensure use of the definition in the `HEMCO_Config.rc.*_metfields` files
- Removed `State_Chm%Spc_Units`

## [14.2.2] - 2023-10-23
### Changed
- Updated sample restart files for fullchem and TransportTracers simulations to files saved out from the 14.2.0 1-year benchmarks

## [14.2.1] - 2023-10-10
### Added
- Script `test/difference/diffTest.sh`, checks 2 different integration tests for differences
- Added GCHP environment file and export/unset env variables in run script for NASA Pleiades cluster
`SatDiagnEdge` collection to all GEOS-Chem Classic `HISTORY.rc` templates
- Added new GCHP config file ESMF.rc for configuring ESMF logging
- Added several new run directory files for use with GEOS-Chem in GEOS
- GCClassic integration tests now display proper commit info in `results.compile.log`
- Stopped OCEAN_CONC from needlessly being pushed through vertical regridding for Hg simulations
- Added warning in GCHP HISTORY.rc about outputting area-dependent variables on custom grids
- Added option to use a single advected species in the carbon simulation
- Added option to perturb CH4 boundary conditions in CH4 simulation
- Added option to perturb OH in CH4 simulation using scale factor in HEMCO_Config.rc

### Changed
- Update `DiagnFreq` in GCClassic integration tests to ensure HEMCO diagnostic output
- Rename restart files in GCHP integration tests (as we do in non-test runs)
- Request 6 hours of execution time for GEOS-Chem Classic integration tests
- Invert directory structure where integration and parallel test scripts are stored
- Error check to stop run if any `MW_g` values are undefined
- Explicitly define tagCH4 simulations in `Input_Opt` rather than basing off of number of advected species
- The `fullchem` mechanism must now be built with KPP 3.0.0 or later
- Changed the AEIC 2019 monthly climatology specification format in ExtData.rc to match standard convention for climatology
- Changed default ESMF logging in GCHP to be ESMF_LOGKIND_NONE (no log)
- NetCDF utilities in `NcdfUtil` folder now use the netCDF-F90 API
- GEOS-only updates for running GEOS-Chem in GEOS
- Boundary conditions for nested-grid simulations are now imposed at every time step instead of 3-hourly
- Update `GeosCore/carbon_gases_mod.F90` for consistency with config file updates in PR #1916
- Update MPI usage in CESM-only code to match new conventions in CAM
- Updated GEPA inventory to GHGI v2 for CH4 and carbon simulations
- Updated integration tests scripts to run on the WashU Compute1 cluster

### Fixed
- Add missing mol wt for HgBrO in `run/shared/species_database_hg.yml`
- Moved the `EDGAR REF_TRF CH4` emissions to the Oil emissions category so it is superseded by GFEIv2 for carbon simulations.
- Prevent `State_Diag%SatDiagnCount` from not being allocated
- For satellite diagnostics, do not test for `id_OH` if OH is not a species
- Fixed parallelization in Luo wetdep simulations caused by uninitialized variable
- Fixed parallelization for Hg0 species in `GeosCore/drydep_mod.F90`
- Fixed incorrect time-slice when reading nested-grid boundary conditions
- Fixed initialization of advected species missing in GCHP restart file
- Fixed comments in `GeosUtil/unitconv_mod.F90` to reflect code implementation
- Fixed compilation issues for `KPP/custom`; updated equations in `custom.eqn`
- Prevent users from creating GCClassic rundirs at 0.25 x 0.3125 resolution for MERRA-2 met
- Added fix to set `RUNDIR_GRID_HALF_POLAR` option for global grids at 0.25x0.3125 or 0.5x0.625 resolutions
- Moved `OCEAN_MASK` out of `ExtData.rc.TransportTracers` and into the
  meteorology template files
- Update `ExtData.rc.CO2` to get meteorology entries from template files
- Added fix for CH4 analytical inversions to convert the state vector value read from file to the nearest integer before comparing to the `Input_Opt%StateVectorElement` read from geoschem_config.yml

### Removed
- Remove references to the obsolete tagged Hg simulation

## [14.2.0] - 2023-10-05
### Added
- Added a printout of GEOS-Chem species and indices
- Added `NcdfUtil/README.md` file directing users to look for netCDF utility scripts at https://github.com/geoschem/netcdf-scripts
- Restored sink reactions for HOI, IONO, IONO2 (fullchem, custom mechanisms)
- S(IV) + HOBr and S(IV) + HOCl reactions to `KPP/fullchem/fullchem.eqn`
- Added setting in GCHP setCommonRunSettings.sh to require species in restarts
- Added setting in GCHP HISTORY.rc to control whether output can be overwritten
- Activated nitrate photolysis
- Added `LightingClimatology` option to HEMCO_Config.rc
- Added run configuration files for WRF-GC
- Added new files `photolysis_mod.F90`, `phot_container_mod.F90`, and `fjx_interface_mod.F90`
- Added photolysis toggle in `geoschem_config.yml` and `Input_Opt` variable Do_Photolysis
- Added speed of light and Planck's constant to PhysConstants module
- Added `GFED4_CLIMATOLOGY` option to HEMCO_Config.rc
- Added CH4 emissions from hydroelectric reservoirs to CH4, Carbon, and tagCH4 simulations
- Added RxnConst diagnostic for archiving reaction rate constants
- Added GCHP run-time option in GCHP.rc to choose dry or total pressure in GCHP advection
- Added GCHP run-time option in GCHP.rc to correct native mass fluxes for humidity
- Added new tracer_mod.F90 containing subroutines for applying sources and sinks for the TransportTracer simulation
- Added new species to the TransportTracer simulation: aoa (replaces CLOCK), aoa_bl, aoa_nh, st80_25
- Added GEOS-IT and GEOSIT as allowable meteorology source options in geoschem_config.yml

### Changed
- Most printout has been converted to debug printout (toggled by `debug_printout: true` in `geoschem_config.yml`
- `HEMCO_Config.rc` template files now use `Verbose: true` to toggle debug printout
- Turn on sea salt debromination via switches in `HEMCO_config.rc`
- If KPP integration fails, reset to prior concentrations and set `RSTATE(3) = 0` before retrying
- Suppress integration errors after 20 errors have been printed to stdout
- Simplified and added comments for bimolecular reactions in clouds in function CloudHet2R
- `HEMCO_Config.rc` and `ExtData.rc` templates now point `HEMCO/GFED4/v2023-03`
- Updated GCHP carbon simulation Global Cl and P(CO) inputs to use 14.0.0 files
- Write GCHP restart files directory to Restarts subdirectory
- Rename GCHP mid-run checkpoint files to standard GEOS-Chem restart format
- Rules for species in restarts files are now the same in GCHP as in GC-Classic
- Moved parts of `CMN_FJX_Mod.F90` not used in original Fast-JX to new container State_Chm%Phot
- Restructured photolysis to create generic photolysis module, interface with Fast-JX, and module containing original Fast-JX analogous to Cloud-J
- Moved UVFlux diagnostics out of JValues collection and into new collection called UVFlux
- In the user registration process:
  - Now ask for both first and last names of the user
  - Now state that user registration is needed for GEOS-Chem support
- Updated `HEMCO_Config.rc` templates to read HEMCO restarts from the `Restarts` rundir folder
- In fullchem simulations, set CO2 to 421 ppm (avg global conc in 2022) everywhere
- Updated CH4 simulation to use CH4 loss frequencies from GCClassic 14.0.0 10-year benchmarks instead of GMI
- Updated CH4 global anthropogenic emission inventory from EDGARv6 to EDGARv7
- Updated `AUTHORS.txt` for version 14.2.0
- Updated links in `README.md` to point to `http://geos-chem.org`
- Changed GCHP default settings to use dry pressure rather than total pressure in advection and correct native mass fluxes for humidity
- Updated partitions requested in Harvard run script examples
- Change RTOL value from 0.5e-3 back to 0.5e-2 to address model slowdown
- Allow the use of OFFLINE_SEASALT for seasalt alkalinity, Cl, and Br in GEOS-Chem within CESM
- Renamed TransportTracer species for consistency with GMAO's TR_GridComp
- See `KPP/fullchem/CHANGELOG_fullchem.md` for fullchem-mechanism
  changes
- Update template `HEMCO_Config.rc.carbon` files to allow running the carbon simulation with only a single species.

### Fixed
- Fixed typo in `GCClassic/createRunDir.sh` preventing benchmark run script from being copied to the run directory
- Fixed divide by zero bug in sulfur chemistry introduced in 14.1.0
- Updates for 0.5 x 0.625 degree GCClassic integration & parallel tests
  - Use `CYS` in `HEMCO_Config.rc` so that missing species in `GC_BCs` will not stop simulations
  - Tests now run for 20 model minutes instead of an hour
- Fixed divide by zero bug in sulfur chemistry introduced in 14.1.0
- Fixed GCHP `HISTORY.rc` issue preventing running with over 3000 cores
- Fixed GCHP `ExtData.rc` error in tagged ozone simulation
- Fixed GCHP `HISTORY.rc` issue preventing diagnostic file overwrite
- Update GCHP interactive run script to fix error handling silent bugs
- Rewrote subroutine calls in `carbon_mod.F90` and `seasalt_mod.F90` to prevent array temporaries.
- Prevent repeated printing of KPP integrate errors to the stdout stream.
- Fixed selection of troposphere-stratosphere boundary in `global_ch4_mod.F90`
- Removed operator splitting in CH4 simulation that was biasing diagnostics
- Fixed GCHP start and elapsed times in time_mod.F90 to use cap_restart value
- Disabled SpeciesConcMND output for benchmark simulations
- Exit `Init_Photolysis` before calling `Calc_AOD` when doing dry-run simulations
- Make sure `State_Het%f_Alk_SSA` and `State_Het%f_Alk_SSC` are in the range 0..1
- Restore seasalt alkalinity to heterogeneous acid-catalyzed reactions of halogens on seasalt aerosols

### Removed
- `Warnings: 1` is now removed from `HEMCO_Config.rc.*` template files
- Removed the `NcdfUtil/perl` folder
- Removed `X-HRS` output from log file
- IONO2 recycling (fullchem, custom mechanisms)
- Deleted unused file set_prof_o3.F90

### Fixed
- Fixed entries for CH4 emissions in `HEMCO_Config.rc.carbon`

## [14.1.2] - 2023-10-20
### Added
- CESM-only update: Added option for correctConvUTLS for correcting buildup of soluble tracers in the UT/LS to match CAM-chem behavior

### Changed
- CESM-only update: extend existing KppError, KppStop to CESM for model stability
- CESM-only update: Removed mpi_bcast in ucx_mod NOXCOEFF_INIT to be handled at coupler level to support spectral-element dynamical core

## [14.1.1] - 2023-03-03
### Added
  - New integration test functions in `test/GCClassic/integration` and `test/GCHP/integration`
  - New parallelization test functions in `test/GCClassic/parallel`
  - Added `README.md` files for integration and parallelization tests in the `test` folder structure
  - Added GCHP integration test for the tagO3 simulation
  - Added GCHP and GCClassic integration tests for the carbon simulation
  - Integration and parallelization test folders have been separated into subdirectories to minimize clutter.
  - GEOS-only updates
  - Add `about` to GitHub issue templates (ensures they will be displayed)
  - Added `.github/ISSUE_TEMPLATE/config.yml` file w/ Github issue options

### Changed
  - GCClassic integration tests now use a single set of scripts
  - GCHP integration tests now use a single set of scripts
  - Integration test run directories are created with the default names assigned by `createRunDir.sh`
  - Several bash functions in `test/shared/commonFunctionsForTests.sh` have been combined so that they will work for both GCClassic and GCHP integration tests
  - `./cleanRunDir.sh` functions now take an argument for non-interactive execution (facilitates integration & parallelization tests)
  - Moved several module variables from `GeosCore/ucx_mod.F90` to `Headers/state_chm_mod.F90`.  This facilitates using GEOS-Chem in CESM.
  - Time cycle flags EFYO are changed to CYS for all GCClassic integration/parallel tests, and for GCClassic fullchem_benchmarksimulations.
- Ask users for the name of their research institution at registration
- Ask users for the name of their PI at registration
- Do not compile GCHP for tagO3 integration tests; use the default build instead
- Moved GC-Classic sample run scripts to operational_examples/harvard_cannon
- The GitHub PR template is now named `./github/PULL_REQUEST_TEMPLATE.md`

### Fixed
- Fixed bug in where writing species metadata yaml file write was always attempted
- Prevent a warning from being generated when compiling `gckpp_Jacobian.F90`
- Fixed a bug in routine GET_IJ where X and Y were swapped in an IF comparison.
- Fixed bug in GFAS pFe by applying work-around in config files

### Removed
- Removed `intTest*_slurm.sh`, `intTest_*lsf.sh`, and `intTest*_interactive.sh` integration test scripts
- Removed State_Met%LWI and input meteorology LWI from carbon simulation run config files
- Removed function `CLEANUP_UCX`; deallocations are now done in `state_chm_mod.F90`

## [14.1.0] - 2023-02-01
### Added
- Added dry deposition updates to Hg0 from Feinberg22 ESPI publication + AMAP emissions
- Added MO2 + NO3 = NO2 + CH2O + HO2 reaction
- Added capability to write species metadata to YAML file
- Added satellite diagnostic (SatDiagn) collection, to archive several fields within a user-defined local-time interval. CAVEAT: For now, only one local-time interval is permitted.
- Added adaptive solver (`rosenbrock_autoreduce`) option for fullchem mechanism
- Added entries for BALD, BENZP, BZCO3H, NPHEN to JValues collection in HISTORY.rc for GCHP
- Added GCHP run script and environment files for MIT clusters Hex and Svante
- Added operational GCHP and GCClassic environment and run scripts for the University of York cluster, Viking
- Added tagO3 run directory for GCHP
- Added upwards mass flux diagnostic to GCHP History collection LevelEdgeDiags
- Added timestep menu to GCHP `geoschem_config.yml` template files
- Added HTAPv3 inventory as a global emissions option (off by default)
- Added carbon simulation and KPP mechanism for CO-CO2-CH4-OCS
- Added GCHP run script and environment file for UCI Australia cluster Gadi
- Added GFAS entries in GCHP config file ExtData.rc

### Changed
- Moved in-module variables in global_ch4_mod.F90 to State_Chm
- Moved in-module variables in hco_interface_gc_mod.F90 to State_Met and State_Chm
- Modified SpeciesConc diagnostic to include option to output units in v/v or molec/cm3
- Rebuilt fullchem and Hg mechanisms with KPP 3.0.0
- Changed HEMCO timecycle flag for QFED and offline emissions from EF to EFY
- Updated the time refresh settings for `O3_PROD` and `O3_LOSS` in `ExtData.rc.tagO3` to read data on the first of each month.

### Fixed
- Fixed sign of Arrhenius "A" coefficient in reaction ETO = HO2 + 2CH2O
- Fixed products in HOBr + SO2 and HOCl + SO2 reactions
- Changed MW_g value of CH4 from 16.05 to 16.04
- Added "WD_CoarseAer:true" for SO4s and NITs in species_database.yml
- Fixed bug in computing State_Met surface type logicals (IsLand, IsWater, etc)
- Fixed bug where State_Met%FRSNO (fraction snow) was all zeros in GCHP
- Fixed HCFC141b and HCFC142b names in GCHP HISTORY.rc
- Fixed list of complex SOA species checked in input_mod.F90
- Now use a string array for reading the list of ObsPack diagnostic species (in `GeosCore/input_mod.F90`)
- Fixed bug in logic that caused restart files not to be linked to the Restarts/ folder of the GCHP tagO3 run directory
- Fixed timestamp for GCClassic History diagnostic so time-averaged collections match the reference time
- Fixed double-titration of seasalt alkalinity
- Fixed bug in GFAS pFe by applying work-around in config files

### Removed
- Removed LRED_JNO2 and AERO_HG2_PARTITON switches from HEMCO_Config.rc (and related code)

## [14.0.2] - 2022-11-29
### Fixed
- Added fix for writing dry-run header to log file
- Updated KPP diagnostics archive flags
- Rewrote code to avoid memory leaks (identified by the code sanitizer)
- Updated EDGAR v6 CH4 emission files to correct timestamp issue
- Updated CH4 Lakes emission files to correct time unit issue
- Added fix for CH4_RICE emissions from EDGAR v6
- Fixed indentation error in the `legacy_bpch` section of `geoschem_config.yml` template files
- Removed "dry air" from the metadata of fields `State_Met%AIRVOL` and `State_Met%BXHEIGHT`
- Applied fixes for CESM runs: Turned off sea salt emissions; Modified time cycle flag for YUAN_MODIS_LAI

### Changed
- Updated CESM HISTORY.rc to work with new CESM-GC diagnostics interface
- Updated sample fullchem restart files copied to run directories to 14.0.0 10-year benchmark output


### Changed
- Use met-field surface type fractions instead of input land-water-ice (LWI) index

### Removed
- Removed State_Met%LWI and LWI as a met-field input

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
- Updated sample restart files copied to run directories to 14.0.0 1-year benchmark output

### Removed
- Removed TMPU1, SPHU1, PS1_WET, and PS1_DRY from GC-Classic restart file
- Removed input.geos; replaced with geoschem_config.yml
- Removed HEMCO.log output file; HEMCO log info now sent to main GEOS-Chem log
