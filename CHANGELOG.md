# Changelog

This file documents all notable changes to the GEOS-Chem repository starting in version 14.0.0, including all GEOS-Chem Classic and GCHP run directory updates.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - TBD
### Added
- Added entries for FINNv25 biomass burning emissions to template HEMCO configuration files
- Added comments to `HEMCO_Diagn.rc` template files instructing users on which ExtNr/Cat/Hier to use for online vs. offline biomass burning emissions
- Added subroutine `Print_Species_Global_Mass` to print_mod for use by GC-Classic
- Added log print of species global mass at start of each timestep if verbose is true
- Added print of global mass computed from restart file values if delta pressure present in restart file
- Added the capability for GCHP simulations to use CH4 restarts for Jacobian Tracers
- Added operational run scripts for WashU Compute2
- Added the option for LPJ_MERRA2 wetland CH4 emissions in CH4 and carbon simulations
- Added GC-Classic config file option to read restart file as `REAL*8` via GEOS-Chem rather than HEMCO
- Added new GCHP run-time option in GCHP.rc to print species mass proxy (Species 1 only) to log file from FV3
- Added GEOS-Chem export in GCHP to send restart file (internal state) delta pressures to FV3 for mixing ratio scaling upon start-up
- Added chemistry budget diagnostics to GCHP carbon HISTORY.rc

### Changed
- Replaced comments in template HEMCO configuration files directing users to obsolete wiki documentation with comments directing users to `hemco.readthedocs.io`
- Updated `EmisOCS_Bioburn` to `EmisOCS_BiomassBurn` in both GCHP `HEMCO_Diagn.rc.carbon` and `HISTORY.rc.carbon` template files
- Updated the ESMF version from 8.4.2 to 8.6.1 in sample environment file `gchp.gcc12_openmpi4_cannon_rocky.env`
- Changed call to `Accept_External_Date_Time` to also pass the seconds value, in order to prevent a WRF-GC bug
- Removed convective washout for default scheme but keep it for LUO_WETDEP
- Adapted Luo2023 WetDep for GF convection
- Updated timestep scaling for convective precipitation areal fraction
- Wrapped tests for infinity/NaN in `#ifdef DEBUG` blocks in `DO_GF_CLOUD_CONVECTION`
- Changed optional argument `Update_Mixing_Ratio` in subroutine `Airqnt` to False by default
- Change GC-Classic call to `Airqnt` to only update mixing ratios if advection is turned off
- Updated mass flux and courant number import scaling in GCHP for compatibility with horizontal flux regridding in MAPL 2.59
- Updated operational run script sample for WashU Compute1
- Update GCHP AWS EFA operational run script examples to avoid crashes over large core counts
- Updated GFEIv3 files to correct issue in original version
- Updated `download_data.py` for compatibility with 0.125 x 0.15625 grids plus all pre-defined nested-grids
- Restructured `download_data.py` to avoid several instances of repeated code
- Changed `read_restart_as_real8` from `false` to `true` in `geoschem_config.yml` for GC-Classic benchmark simulations
- Changed the default setting of `read_restart_as_real8` from `false` to `true` in template file `geoschem_config.yml.TransportTracers`
- Disable PARANOX extension when using GEOS-Chem Classic 0.25x0.3125 or 0.125x0.15625 grids
- Commented out met-fields `PEDGEDRY`, `PFICU`, `PFILSAN`, `PFLCU`, and `PFLLSAN` by default in GC-Classic and GCHP carbon HISTORY.rc, and GC-Classic CH4 HISTORY.rc
- Turned on Carbon collection in `HISTORY.rc` for carbon simulations by default

### Fixed
- Restored entries for TMB emissions in `HEMCO_Config.rc.fullchem` template files for GCClassic and GCHP
- Moved `EmisOCS_Total` to the head of the `EmisOCS` diagnostic entries in the GCHP `HISTORY.rc.carbon` template file
- Fixed OM/OC ratio for OCPO in SimpleSOA to be 1.4 instead of 2.1
- Fixed precipitation formation rate unit in Luo2023 convective washout
- Fixed bug where species mass in restart file was not conserved in first timestep if run-time meteorology different from restart file meteorology
- Fixed parallel errors in `convection_mod.F90` by setting `AER = . TRUE.` and `KIN = .TRUE.` before calling `WASHOUT`
- Fixed Hg directional ocean flux diagnostics in the Hg simulation so that they equal net flux
- Fixed error where `//` were not being changed to `/` in `download_data.py`
- Change precision of area import from GCHP advection from `REAL*4` to native `REAL*8`
- Fixed time-range and units for CH4 emission inventories to be consistent with the corresponding netCDF files in ExtData directory for `HEMCO_Config.rc` and `ExtData.rc`
- Updated scaling factor ID at 3000 to avoid conflicts with CEDS_01x01 scaling factor enabled in carbon simulation for IMI analytical inversion
- Fixed typos in `ind_` variable names in `KPP/carbon/carbon_Funcs.F90`
- Fixed typo in GCHP operational run script for Harvard Cannon to properly retrieve the run duration string
- Fixed bug in ObsPack to include instantaneously-sampled data whose
  timestamps are within 1/2 of a model timestep of the end of the day
- Fixed out-of-bounds error in `carbon_gases_mod.F90` that is caused by refernencing `OHdiurnalFac` array when it is not defined

### Removed
- Removed entries for FINN v1.5 biomass burning emissions from template HEMCO configuration files
- Removed `Is_Advected` tags from `run/shared/species_database*.yml` template files
- Removed GCHP initialization of `State_Met` fields `TropLev`, `BxHeight`, and `DELP_DRY` from restart file values since over-written with values of current meteorology
- Removed `OH_PosteriorSF` entry in carbon and CH4 HEMCO_Config.rc since never used
- Removed extraneous division by `TS_EMIS` in routine `Chem_H2O2` (located in `GeosCore/sulfate_mod.F90`)

## [14.6.3] - 2025-07-28
### Added
- Added error check to exclude sampling ObsPack observations located outside of a nested-grid domain
- Added Grell-Freitas convection subroutine for post-GEOS-5.22 (GEOS-IT and GEOS-FP after June 2020)
- Added GEOS-IT simulations to use offline emissions generated with GEOS-IT
- Added meteorology-specific `OFFLINE_EMISSION_DIR` entries in shared directory for future use and backward compatibility
- Added operational run script sample for AWS with EFA-enabled
- Added operational run scripts for Harvard Cannon with Intel VTune commands
- Added sample environment file for Harvard Cannon with GNU 14.2.0 compilers
- Converted `F` in `DO_CONVECTION` from a variable to a pointer, for computational speedup
- Changed OpenMP loop scheduling from `DYNAMIC` to `GUIDED` in routine `DO_CONVECTION`
- Added `Diagn_APM` routine in `GeosCore/hcoi_gc_diagn_mod.F90` to restore HEMCO manual diagnostics for use w/ APM
- Added hidden option to read GC-Classic restart file as real8 locally rather than real4 through HEMCO

### Changed
- Updated logic to include ObsPack observations that span UTC date boundaries
- Assigned ObsPack averaging interval end times (instead of start times) to the `aveEnd` variable in routine `ObsPack_Write_Output`
- Optimized parallel loops in `AIRQNT` routine in `GeosCore/calc_met_mod.F90`
- Optimized parallel loops in `VDIFF` routine in `GeosCore/vdiff_mod.F90`
- Placed error checks for infinity or NaN in `DO_CONVECTION` in `#ifdef DEBUG` preprocessor blocks
- Collapsed several parallel DO loops in `GeosCore/carbon_mod.F90`
- Changed met guidance in run directory creation to remove beta for GEOS-IT, make GCHP mass fluxes beta, and improve GEOS-FP warning
- Changed path to carbon, CH4, CO2 simulation restart files to `ExtData/GEOSCHEM_RESTARTS/v2025-07` in `download_data.yml` and the GCHP `createRunDir.sh` script

### Fixed
- Added missing 3rd element in assigment of `Item%NcChunkSizes` in `History/histitem_mod.F90`
- Removed extra unit conversion to mol/mol on 0th hour boundary conditions in `History/history_mod.F90`
- Reordered code in `aerosol_mod.F90` and `gc_environment_mod.F90` so that aerosol optics file paths will be printed to the dry-run log file
- Fixed wrong mass flux and Courant number import scaling for GCHP runs that read these fields from offline files
- Corrected GCHP carbon HISTORY.rc entries for KPPdiags, RxnRates, and RxnConst collections

### Removed
- Removed `#ifndef TOMAS` block at the start of the parallel loop in `DO_CONVECTION`
- Removed redundant `IF/ELSE` statement in the 2nd parallel loop in routine `AIRQNT`
- Removed redundant `ELSE` blocks in `DO_CONVECTION`
- Removed redundant `units` variable in routine `AIRQNT`

## [14.6.2] - 2025-06-11
### Added
- Added MCHgMAP geogenic emissions (2010-2020) from Dastoor et al. (2025)
- Added MCHgMAP biomass emissions (2010-2020, GFED and FINN) from Dastoor et al. (2025)
- Added NCAR Derecho operational example environment file `geoschem.intel24.env`

### Changed
- Modified `ObsPack/obspack_mod.F90` to use GEOS-Chem surface geopotential height (`PHIS`) for selecting model layer for comparison to obspack.
- Updated NCAR Derecho run script to source the `geoschem.intel24.env` environment file

### Fixed
- Restored unit convertion for boundary conditions from mol/mol to kg/kg dry air

### Fixed
- Assigned ObsPack averaging interval end times (instead of start times) to the `aveEnd` variable in routine `ObsPack_Write_Output`
- Fixed logic in `Obspack_Sample` to include observations whose averaging windows span a UTC date boundary

## [14.6.1] - 2025-05-27
### Added
- Added `#InvCEDSshipALK6`, `#InvCEDS_TMB`, and `#InvCEDSship_TMB` (commented out by default) to `HEMCO_Diagn*` and GCHP `HISTORY.rc.fullchem` files
- Added `EmisOCs*` diagnostics to `run/GCHP/HEMCO_Diagn.rc.templates/HEMCO_Diagn.rc.carbon` (these were missing)
- Added entry for GFAS methanol for ExtData
- Added `RxnConst` and `RxnRates` History collections to `HISTORY.rc.carbon` and `HISTORY.rc.Hg` template files
- Added routine `Hg_UpdateKppDiags` to update the `KppDiags` history diagnostic arrays in `mercury_mod.F90`
- Added `rrtmg_radiative_transfer_model:aod_wavelength_in_nm` to the GCClassic `geoschem_config.yml.aerosol` template file
- Added an error trap in routine `Init_Aerosol` to make sure at least 1 AOD wavelength is selected for fullchem or aerosol-only simulations
- Added operational example run script for GEOS-Chem Classic on NCAR Derecho cluster
- Added options to run GCHP using C720 mass fluxes or derived winds with 0.25x0.3125 processed files for other meteorology

### Changed
- Updated GCHP template files `HEMCO_Diagn.rc.fullchem` and `HISTORY.rc.fullchem` so that the same emission diagnostics are requested in both
- Changed `ALD2_PLANTDECAY` emissions category (for GEOS-Chem in NASA-GEOS ESM only) from 3 to 99 to not conflict with the anthropogenic transport sector
- Abstracted diagnostic code out of `Chem_Carbon_Gases` and into PRIVATE subroutines in `GeosCore/carbon_gases_mod.F90`
- Abstracted diagnostic code out of `ChemMercury` and into PRIVATE subroutines in `GeosCore/mercury_mod.F90`
- Modified logic in `Init_State_Diag` so `KppDiags` diagnostic fields can be registered when using fullchem, Hg, or carbon mechanisms
- Modified logic in `Init_State_Diag` so that `JValues` and `UVFlux` diagnostic fields can be registered when using fullchem or Hg simulations
- Modified `Obspack_Read_Input` routine to look for `middle of the averaging interval` and `midpoint of the averaging interval` in the `time:comment` string
- Wrapped several TOMAS print statements in `IF ( Input_Opt%Verbose )` blocks to avoid excessive printout when using GCHP-TOMAS
- Changed GCHP recommended GEOS-IT options for meteorology from mass fluxes with raw C180 fields to 3hr winds with processed C180 fields.

### Fixed
- Restored the `UVFlux` diagnostic collection to the GCHP `fullchem_alldiags` integration test
- Fixed outdated path for GFED4 daily fraction
- Fixed entries for GEOS-IT preprocessed cubed-sphere wind in GCHP
- Fixed the `KppTime` diagnostic in `Chem_Carbon_Gases`; it was not being updated properly

### Removed
- Removed unused run directory creation files for GCHP grid resolutions c24 and c48
- Fixed index-based access of hydrophobic and hydrophilic carbon species in TOMAS.
- Fixed the `KppTime` diagnostic in `Chem_Carbon_Gases`; it was not being updated properly
- Fixed a bug in GCClassic and GCHP integration test scripts that caused `KppTime` not to be commented out in `fullchem_alldiags` tests
- Removed most entries under the `photolysis` section in `geoschem_config.yml.aerosol` template file, as the aerosol-only simulation doesn't call Cloud-J
- Removed setting `DELP_DRY` to zero in `hco_utilities_gc_mod.F90` when not found in the restart file to avoid negative concentrations

## [14.6.0] - 2025-04-18
### Added
- Added CEDS 0.1 x 0.1 degree emissions (in `HEMCO/CEDS/v2024-06`)
- Added met-field dependent dust tuning factors for GCHP C24 resolution in `setCommonRunSettings.sh`.  Others to be added later.
- Added placeholder values for dust mass tuning factors in `HEMCO_Config.rc.GEOS`
- Added dust scale factors for MERRA-2, GEOS-IT, and GEOS-FP when using USTAR for Dust DEAD extension
- Added utility subroutine `Print_Species_Min_Max_Sum` to `print_mod.F90`
- Added routine `Set_DryDepVel_Diagnostics` to `hco_interface_gc_mod.F90`
- Added dry-run integration tests for selected simulations
- Added `State_Diag%SatDiagnPMid` and `State_Diag%Archive_SatDiagnPMid` to save pressure at level midpoints to the `SatDiagn` collection
- Added option to run GCClassic nested-grid simulations at 0.125x0.15625 resolution using GEOS-FP derived winds fields generated from c720 mass fluxes archived by GMAO
- Added option for South America (SA), Africa (AF), Middle East (ME), Oceania (OC), and Russia (RU) regions to nested-grid simulations in GCClassic's createRunDir.sh
- Added updates for compatibility with the Beijing Climate Centre Earth System Model
- Added `run/shared/rtd_species_by_simulation.py` script to generate tables of species for each simulation for ReadtheDocs

### Changed
- Updated default CEDS from CEDSv2 (0.5 deg x 0.5 de) to new CEDS (0.1 deg x 0.1 deg)
- Added the `KPP_INTEGRATOR_AUTOREDUCE` C-preprocessor switch integrator-specific handling
- Added code to `KPP/*/CMakeLists.txt` to read the integrator name from the `*.kpp` file
- Replaced `GOTO` statements with `IF/THEN/ELSE` blocks in `GeosCore/drydep_mod.F90`
- Changed several diagnostic subroutines to expect species concentrations in mol/mol rather than kg/kg
- Added precision when registering `State_Met` and `State_Chm` arrays to change output file precision to match precision in the model
- Changed GEOS-Chem Classic restart file precision of species concentrations (`State_Chm%SpcRestart`) from `REAL*4` to `REAL*8` to match precision in the model
- Moved GEOS-Chem Classic retrieval of restart variable DELPDRY from HEMCO to `GC_Get_Restart` for consistency with handling of all other restart variables
- Moved dry dep velocity diagnostic outputs for sea flux and satellite diagnostic species into the `DryDep` collection
- Moved computation of the `DryDepVelForAlt1` diagnostic into routine `Set_DryDepVel_Diagnostics` 
- Updated `run/shared/download_data.yml` to use `--no-sign-request` for S3 downloads via anonymous login
- Changed `KPP/CMakeLists.txt` to not call `add_directory(standalone)` unless we have configured with `-DKPPSA=y`
- Moved Cloud-J and Fast-JX input directories to Cloud-J and new Fast-JX menus respectively in `geoschem_config.yml`
- Updated photolysis and aerosol optics input directories to use new mineral dust values in `FJX_scat-aer.dat` and `dust.dat` based on spheroidal shapes
- Set `State_Diag%Archive_SatDiagn` to true if `State_Diag%Archive_SatDiagnPMID` is true
- Updated `RxnRates` and `RxnConst` diagnostic fields to use 4-digit reaction numbers.
- Rebuilt `fullchem`, `Hg`, `carbon` chemical mechanisms with KPP 3.2.0
- Changed the minimum KPP version to 3.2.0
- Disabled the `KppTime` diagnostic output in the `fullchem_alldiags` integration tests; this will vary from run to run causing difference tests to fail
- Updated the `KPP-Standalone` for compatibility with KPP 3.2.0 and to write the proper number of header lines to skip before data begins
- Set `use_archived_PCO_from_CH4` and `use_archived_PCO2_from_CO2` to true by default for carbon simulations
- Updated CH4 global oil, gas, and coal emissions from GFEIv2 to GFEIv3
- Changed GCHPctmEnv and DYNAMICS diagnostic names in GCHP to include suffix '_R4'

### Fixed
- Fixed PDOWN definition to lower rather than upper edge
- Moved where prescribed CH4 is applied in GEOS-Chem Classic to after emissions application so that updated PBL heights are used
- Moved species concentration unit conversions between mol/mol and kg/kg to start and end of every timestep in GEOS-Chem Classic to remove differences introduced when reading and writing restart files
- Fixed bug in restart file entry for `ORVCSESQ` in GEOS-Chem Classic fullchem HEMCO_Config.rc that resulted in initializing to all zeros
- Fixed parallelization issue when computing `State_Chm%DryDepNitrogren` used in HEMCO soil NOx extension
- Fixed bugs in column mass array affecting budget diagnostics for fixed level and PBL
- Updated `SatDiagnColEmis` and `SatDiagnSurfFlux` arrays in `hco_interface_gc_mod.F90`, with `(I,J,S)` instead of `(:,:,S)`
- Fixed incorrect description metadata for `FluxHg0FromAirToOcean` and `FluxHg0FromOceanToAir` diagnostics
- Placed call to `Convert_Spc_Units` in `main.F90` within an `IF ( notDryrun )` block to avoid executing unit conversions in dryrun simulations
- Modified CH4 reservoir timestamps in HEMCO_Config.rc to use months 1-12 to ensure HEMCO recalculates those fields monthly and properly applies the seasonal mask
- Fixed path error `download_data.py` when downloading from `geoschem+http` or `nested+http` portals
- Retrieve UV flux arrays from Cloud-J used to set UV flux diagnostics
- Fixed issue in `download_data.py` that was adding an extra `ExtData` to file paths
- Restored `UVFlux` diagnostic output in the `fullchem_alldiags` integration test
- Restored convection and ConvertBox unit conversion parallelization to how it was prior to 14.4.0 to fix slowness in TOMAS simulations
- Modified the carbon mechanism in KPP to separate tropospheric CH4 loss by OH from CO production by CH4 to remove dependency of CH4 and CO on each other and eliminate differences between CH4/tagCO simulations and the carbon simulation
- Renamed several dummy species in the carbon mechanism for clarity
- Fixed precision calculations within `co2_mod.F90` and `tagged_co_mod.F90` to eliminate differences with the carbon simulation
- Fixed simulation date information printed by metrics.py for GCHP

### Removed
- Removed `CEDSv2`, `CEDS_GBDMAPS`, `CEDS_GBDMAPSbyFuelType` emissions entries from HEMCO and ExtData template files
- Removed re-evaporation requirement for washout
- Removed unused level argument passed to `SOIL_DRYDEP` and `SOIL_WETDEP`
- Removed Fast-JX input directory from geoschem_config.yml files except for Hg simulation
- Removed `History` attribute from ObsPack output netCDF files; the date info was causing difference tests to fail
- Removed unused diagnostics: `Tomas_H2SO4`, `Tomas_COAG`, `Tomas_NUCL`, `Tomas_AQOX`, `Tomas_MNFIX`, `Tomas_SOA`
- Removed diurnal cycle factor applied to OH in `KPP/carbon/carbon_Funcs.F90` to eliminate differences between CH4 and carbon simulations.
- Removed diurbal cycle factor applied to OH in tagCO simulation for consistency with other carbon species
- Removed unused functions from `carbon_get_CO2fromOH_flux` and `carbon_get_FixedOH_Flux` from `KPP/carbon/carbon_Funcs.F90`

## [14.5.3] - 2025-03-04
### Changed
- Changed CESM `HEMCO_Config.rc` to read 3D AEIC emissions every timestep to avoid differences upon restart

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
- Added Australian Hg emissions for 2000-2019 from MacFarlane et. al. [2022], plus corresponding mask file
- Added comments in GEOS-Chem Classic `HISTORY.rc` template files advising users not to change the `BoundaryConditions.frequency` setting
- Added `.zenodo.json` for auto-DOI generation upon version releases

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
- Reverted CH4 livestock emissions to EDGAR v7 to avoid hotspots and to apply seasonality

### Removed
- Removed duplicate `WD_RetFactor` tag for HgClHO2 in `species_database.yml`
- Removed error messages in HEMCO interface pointing users to HEMCO log
- Removed unused RUNDIR settings for GCHP pressure units and scaling

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
