add_library(History STATIC
	histcontainer_mod.F90
	histitem_mod.F90
	history_mod.F90
	history_netcdf_mod.F90
	history_util_mod.F90
	metahistcontainer_mod.F90
	metahistitem_mod.F90
)
target_link_libraries(History
	PUBLIC KPP
)

