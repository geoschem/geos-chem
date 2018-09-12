#[[
Uses variables:
    MECH
]]

if("${MECH}" STREQUAL "Standard")
    add_subdirectory(Standard)    
else()
    message(FATAL_ERROR "Build system for MECH, ${MECH}, is not yet supported.")
endif()