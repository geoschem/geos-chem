# Headers/CMakeLists.txt

#----------------------------------------------------------------------------
# Define libHeaders.a
#----------------------------------------------------------------------------
add_library(Headers
  STATIC EXCLUDE_FROM_ALL
  charpak_mod.F90
  CMN_DIAG_mod.F90
  CMN_FJX_MOD.F90
  CMN_O3_mod.F90
  CMN_SIZE_mod.F90
  diaglist_mod.F90
  dictionary_m.F90
  taggeddiaglist_mod.F90
  errcode_mod.F90
  input_opt_mod.F90
  inquireMod.F90
  physconstants.F90
  precision_mod.F90
  qfyaml_mod.F90
  registry_mod.F90
  registry_params_mod.F90
  roundoff_mod.F90
  species_database_mod.F90
  species_mod.F90
  state_chm_mod.F90
  state_diag_mod.F90
  state_grid_mod.F90
  state_met_mod.F90
)

# Define dependencies for libHeaders.a
target_link_libraries(Headers PUBLIC KPP_FirstPass)

