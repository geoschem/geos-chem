if("${MECH}" STREQUAL fullchem)
	add_subdirectory(fullchem)
endif()
if("${MECH}" STREQUAL custom)
	add_subdirectory(custom)
endif()

