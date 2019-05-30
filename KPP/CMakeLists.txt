# Add the subdirectory with desired chemistry mechanism
if("${MECH}" STREQUAL "Standard")
    add_subdirectory(Standard)
elseif("${MECH}" STREQUAL "Tropchem")
    add_subdirectory(Tropchem)
elseif("${MECH}" STREQUAL "SOA_SVPOA")
    add_subdirectory(SOA_SVPOA)
else()
    message(FATAL_ERROR "KPP build target for MECH=${MECH} is not yet implemented")
endif()

