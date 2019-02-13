# Add the correct subdirectory based on ${MECH}
if("${MECH}" STREQUAL "Standard" OR "${MECH}" STREQUAL "benchmark")
    add_subdirectory(Standard)    
elseif("${MECH}" STREQUAL "Tropchem")
    add_subdirectory(Tropchem)
else()
    message(FATAL_ERROR "KPP build target for MECH=${MECH} is not yet implemented")
endif()

