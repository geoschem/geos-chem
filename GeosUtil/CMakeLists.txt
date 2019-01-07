add_library(GeosUtil STATIC
	henry_mod.F
	ifort_errmsg.F
	julday_mod.F
)
target_link_libraries(GeosUtil
	PUBLIC Headers
)

add_library(GeosUtil2 STATIC
	bpch2_mod.F
	error_mod.F
	file_mod.F
	gc_grid_mod.F90
	geos_timers_mod.F90
	grid_registry_mod.F90
	hdf_mod.F
	pressure_mod.F
	regrid_a2a_mod.F90
	roundoff_mod.F90
	time_mod.F
	transfer_mod.F
	unitconv_mod.F90
)
target_link_libraries(GeosUtil2
	PUBLIC NcdfUtil
)
