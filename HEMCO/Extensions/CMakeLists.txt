add_library(HCOX STATIC EXCLUDE_FROM_ALL
	drydep_toolbox_mod.F90
	hcox_ch4wetland_mod.F90
	hcox_custom_mod.F90
	hcox_driver_mod.F90
	hcox_dustdead_mod.F
	hcox_dustginoux_mod.F90
	hcox_finn_include.FINNv16
	hcox_finn_mod.F90
	hcox_finn_mod.FINNv16
	hcox_gc_POPs_mod.F90
	hcox_gc_RnPbBe_mod.F90
	hcox_gfed_mod.F90
	hcox_iodine_mod.F90
	hcox_lightnox_mod.F90
	hcox_megan_mod.F90
	hcox_paranox_mod.F90
	hcox_seaflux_mod.F90
	hcox_seasalt_mod.F90
	hcox_soilnox_mod.F90
	hcox_state_mod.F90
	hcox_template_mod.F90x
	$<$<BOOL:${TOMAS}>:hcox_tomas_jeagle_mod.F90 hcox_tomas_dustdead_mod.F>
	hcox_tools_mod.F90
	hcox_volcano_mod.F90
	ocean_toolbox_mod.F90
)
target_link_libraries(HCOX
	PUBLIC HCO
)

