add_library(Headers STATIC EXCLUDE_FROM_ALL
	charpak_mod.F90
	CMN_DIAG_mod.F
	CMN_FJX_MOD.F
	CMN_O3_mod.F
	CMN_SIZE_mod.F
	diaglist_mod.F90
	errcode_mod.F90
	input_opt_mod.F90
	inquireMod.F90
	physconstants.F
	precision_mod.F
	registry_mod.F90
	registry_params_mod.F90
	species_database_mod.F90
	species_mod.F90
	state_chm_mod.F90
	state_diag_mod.F90
	state_grid_mod.F90
	state_met_mod.F90
)
target_link_libraries(Headers
	PUBLIC
	$<$<STREQUAL:${MECH},SOA_SVPOA>:KPPFirstPass_SOA_SVPOA>
	$<$<STREQUAL:${MECH},Standard>:KPPFirstPass_Standard>
	$<$<STREQUAL:${MECH},Tropchem>:KPPFirstPass_Tropchem>
)

