add_library(APM STATIC EXCLUDE_FROM_ALL
	apm_albd_mod.F
	apm_ATHN_mod.F
	apm_coag_mod.F90
	apm_grow_mod.F	
	apm_icen_mod.F90
	apm_init_mod.F
	apm_mixactivate.F90
	apm_nucl_mod.F
	apm_opti_mod.F
	apm_phys_mod.F
	apm_radf_mod.F
	apm_rrtmg_mods.F90
	apm_rrtmg_src.F90
	apm_rrtmg_sw.F90
	apm_TIMN_mod.F
	module_data_mosaic_asect.F90
	module_data_mosaic_other.F90
	module_data_mosaic_therm.F90
	module_mosaic_therm.F90
)
target_link_libraries(APM
	PUBLIC GeosUtil
)

