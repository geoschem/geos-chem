add_library(HCOI STATIC
	hcoi_esmf_mod.F90
	hcoi_standalone_mod.F90
)
target_link_libraries(HCOI
	PUBLIC HCOX
)

if("${HCOSA}")
	add_executable(hemco_standalone
		hemco_standalone.F90
	)
	target_link_libraries(hemco_standalone
		PUBLIC HCOI
	)
	install(TARGETS hemco_standalone
		RUNTIME DESTINATION ${RUNDIR}
	)
endif()
