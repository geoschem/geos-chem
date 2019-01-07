add_library(HCOX STATIC
	drydep_toolbox_mod.F90
	hcox_aerocom_mod.F90
	hcox_ch4wetland_mod.F90
	hcox_custom_mod.F90
	hcox_driver_mod.F90
	hcox_dustdead_mod.F
	hcox_dustginoux_mod.F90
	hcox_finn_include.H
	hcox_finn_mod.F90
	hcox_gc_POPs_mod.F90
	hcox_gc_RnPbBe_mod.F90
	hcox_gfed_include_gfed3.H
	hcox_gfed_include_gfed4.H
	hcox_gfed_mod.F90
	hcox_iodine_mod.F90
	hcox_lightnox_mod.F90
	hcox_megan_mod.F
	hcox_paranox_mod.F90
	hcox_seaflux_mod.F90
	hcox_seasalt_mod.F90
	hcox_soilnox_mod.F90
	hcox_state_mod.F90
	hcox_tools_mod.F90
	ocean_toolbox_mod.F90

	# TODO: can these be omitted?
	# hcox_tomas_dustdead_mod.F
	# hcox_tomas_jeagle_mod.F90
)
target_link_libraries(HCOX
	PUBLIC HCO
)
