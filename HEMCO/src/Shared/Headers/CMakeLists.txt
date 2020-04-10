
add_library(HeadersHco STATIC EXCLUDE_FROM_ALL
	charpak_mod.F90
	inquireMod.F90
	physconstants.F
	precision_mod.F
)
target_link_libraries(HeadersHco
	PUBLIC HEMCOBuildProperties
)
target_include_directories(HeadersHco
	INTERFACE ${HEMCO_BINARY_DIR}/mod
)
