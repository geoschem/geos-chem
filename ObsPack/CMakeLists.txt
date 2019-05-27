add_library(ObsPack
    obspack_mod.F90
)
target_link_libraries(ObsPack
	PUBLIC History
)