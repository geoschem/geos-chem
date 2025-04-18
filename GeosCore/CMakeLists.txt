# GeosCore/CMakeLists.txt

# Add a compilation flag to force REAL to REAL*8 if necessary
# only for file tpcore_window_mod.F90
if(${USE_REAL8})
  if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    set_source_files_properties(tpcore_window_mod.F90
      PROPERTIES COMPILE_FLAGS "-r8")
  elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    set_source_files_properties(tpcore_window_mod.F90
      PROPERTIES COMPILE_FLAGS "-fdefault-real-8 -fdefault-double-8")
  endif()
endif()

#-----------------------------------------------------------------------------
# Define libGeosCore.a
#-----------------------------------------------------------------------------
add_library(GeosCore
  STATIC EXCLUDE_FROM_ALL
  aero_drydep.F90
  aerosol_mod.F90
  aerosol_thermodynamics_mod.F90
  airs_ch4_mod.F90
  calc_met_mod.F90
  carbon_mod.F90
  carbon_gases_mod.F90
  chemistry_mod.F90
  cldice_HBrHOBr_rxn.F90
  cldj_interface_mod.F90
  cleanup.F90
  co2_mod.F90
  convection_mod.F90
  depo_mercury_mod.F90
  diagnostics_mod.F90
  drydep_mod.F90
  dust_mod.F90
  emissions_mod.F90
  fast_jx_mod.F90
  fjx_interface_mod.F90
  fullchem_mod.F90
  gc_environment_mod.F90
  get_ndep_mod.F90
  global_br_mod.F90
  global_ch4_mod.F90
  gosat_ch4_mod.F90
  hcoi_gc_diagn_mod.F90
  hco_state_gc_mod.F90
  hco_interface_gc_mod.F90
  hco_utilities_gc_mod.F90
  input_mod.F90
  kppsa_interface_mod.F90
  land_mercury_mod.F90
  linear_chem_mod.F90
  linoz_mod.F90
  mapping_mod.F90
  mercury_mod.F90
  mixing_mod.F90
  modis_lai_mod.F90
  ocean_mercury_mod.F90
  olson_landmap_mod.F90
  pbl_mix_mod.F90
  photolysis_mod.F90
  pjc_pfix_mod.F90
  pjc_pfix_window_mod.F90
  planeflight_mod.F90
  pops_mod.F90
  RnPbBe_mod.F90
  rpmares_mod.F90
  seasalt_mod.F90
  set_global_ch4_mod.F90
  sfcvmr_mod.F90
  sulfate_mod.F90
  tagged_co_mod.F90
  tagged_o3_mod.F90
  tccon_ch4_mod.F90
  toms_mod.F90
  tpcore_fvdas_mod.F90
  tpcore_window_mod.F90
  tracer_mod.F90
  transport_mod.F90
  ucx_mod.F90
  uvalbedo_mod.F90
  vdiff_mod.F90
  wetscav_mod.F90
  YuIMN_Code.F90

  # Files only included for special cases
  $<$<BOOL:${MODEL_CLASSIC}>:flexgrid_read_mod.F90 get_met_mod.F90 set_boundary_conditions_mod.F90>
  $<$<BOOL:${RRTMG}>:rrtmg_rad_transfer_mod.F90>
  $<$<BOOL:${TOMAS}>:tomas_mod.F90>
  $<$<BOOL:${APM}>:apm_driv_mod.F90>
  # exchange_mod : TODO
)

# Fixes #64: If GCC>=8.0.0 then ocean_mercury_mod.F90 
# should be compiled with -O1
if(("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU") AND NOT (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "8.0.0") AND NOT("${CMAKE_BUILD_TYPE}" STREQUAL "Debug"))
  set_source_files_properties(ocean_mercury_mod.F90
    PROPERTIES COMPILE_FLAGS -O1
  )
endif()

# Define dependencies for libGeosCore.a
target_link_libraries(GeosCore PUBLIC
  ObsPack 
  History
  KPP
  HETP_core
  $<LINK_ONLY:HCOI_Shared>    # link only to avoid gathering the rest of the TURs
  $<$<BOOL:${APM}>:APM>       # depends on APM iff APM is true
  $<$<BOOL:${RRTMG}>:GeosRad> # depends on GeosRad iff RRTMG is true
  $<$<BOOL:${GTMM}>:Hg>       # depends on Hg iff GTMM is true

)
