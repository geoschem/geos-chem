add_library(ObsPack STATIC EXCLUDE_FROM_ALL
    obspack_mod.F90
)
target_link_libraries(ObsPack
	PUBLIC GeosUtil
)