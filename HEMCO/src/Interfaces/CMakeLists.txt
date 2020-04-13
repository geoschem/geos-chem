add_library(HCOI STATIC EXCLUDE_FROM_ALL
	hcoi_esmf_mod.F90
	hcoi_standalone_mod.F90
)
target_link_libraries(HCOI
	PUBLIC HCOX
)
target_compile_options(HCOI
	PRIVATE
		""
		$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-ffree-line-length-none>
)

if("hemco_standalone" IN_LIST GCCLASSIC_EXE_TARGETS)
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

