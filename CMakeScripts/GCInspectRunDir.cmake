#[[ GCInspectRunDir

This module inspects the run directory and sets the following variables:
    RUNDIR_MET
    RUNDIR_NESTED
    RUNDIR_REGION
    RUNDIR_GRID
    RUNDIR_MECH

This is done by executing the getRunInfo script in the run directory. Note
that the responses of that function, e.g. "geosfp", are modified to the 
expected values of GCBuildFlags, e.g. "GEOS_FP".
]]

message(STATUS "Bootstrapping ${RUNDIR}")

# Make sure that the run directory is defined
if(NOT DEFINED RUNDIR)
    message(FATAL_ERROR "RUNDIR not defined")
endif()

macro(inspect_rundir VAR ID)
    execute_process(COMMAND perl ${RUNDIR}/getRunInfo ${RUNDIR} ${ID}
        OUTPUT_VARIABLE ${VAR}
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
endmacro()

# Inspect MET
inspect_rundir(RUNDIR_MET 0)
if("${RUNDIR_MET}" STREQUAL "geosfp")
    set(RUNDIR_MET "GEOS_FP")
elseif("${RUNDIR_MET}" STREQUAL "merra2")
    set(RUNDIR_MET "MERRA2")
endif()

# Inspect GRID
inspect_rundir(RUNDIR_GRID 1)
if("${RUNDIR_GRID}" STREQUAL "2x25")
    set(RUNDIR_GRID "2x2.5")
elseif("${RUNDIR_GRID}" STREQUAL "05x0625")
    set(RUNDIR_GRID "0.5x0.625")
elseif("${RUNDIR_GRID}" STREQUAL "025x03125")
    set(RUNDIR_GRID "0.25x0.3125")
endif()

# Inspect MECH
inspect_rundir(RUNDIR_MECH 2)
if("${RUNDIR_MECH}" STREQUAL "standard")
    set(RUNDIR_MECH "Standard")
elseif("${RUNDIR_MECH}" STREQUAL "tropchem")
    set(RUNDIR_MECH "Tropchem")
endif()

# Inspect NESTED
inspect_rundir(RUNDIR_REGION 3)
if("${RUNDIR_REGION}" STREQUAL "n")
    set(RUNDIR_NESTED "FALSE")
    unset(RUNDIR_REGION)
else()
    set(RUNDIR_NESTED "TRUE")
    string(TOUPPER "${RUNDIR_REGION}" RUNDIR_REGION)
endif()

