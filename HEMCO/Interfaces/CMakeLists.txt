add_library(HCOI STATIC
	hcoi_esmf_mod.F90
	hcoi_standalone_mod.F90
	hemco_standalone.F90
)
target_link_libraries(HCOI
	PUBLIC HCOX
)
