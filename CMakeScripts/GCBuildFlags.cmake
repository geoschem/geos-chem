#[[ GCBuildFlags.cmake

This module builds the list of compiler options and definitions for 
building GEOS-Chem. Note that compiler options are platform and compiler
dependent, whereas definitions depend on the build configuration (i.e. MET, 
GRID, etc.).

Exported variables:
    IMPL
    GC_DEFINES
    FC_OPTIONS

    RRTMG

]]


# Include logging and setting helpers
include(GCHelpers)

# Guard against multiple includes
if(DEFINED GC_BUILD_FLAGS_INCLUDED)
    return()
else()
    set(GC_BUILD_FLAGS_INCLUDED "TRUE")
endif()

# Get run directory
message(STATUS "Run directory setup:")

set(RUNDIR_PROMPT "<path to run directory>")
set_dynamic_default(RUNDIR "${RUNDIR_PROMPT}"
    LOG RUNDIR_LOG
)
dump_log(RUNDIR_LOG)

message(STATUS "Bootstrapping ${RUNDIR}")

# Check that RUNDIR exists
if("${RUNDIR}" STREQUAL "${RUNDIR_PROMPT}")
    message(FATAL_ERROR "You haven't set RUNDIR!")
elseif(NOT IS_ABSOLUTE "${RUNDIR}")
    set(RUNDIR "${CMAKE_BINARY_DIR}/${RUNDIR}")
endif()

if(NOT EXISTS "${RUNDIR}")
    message(FATAL_ERROR "Invalid RUNDIR. ${RUNDIR} does not exist!")
endif()

# Inspect the RUNDIR's name to make a guess at if we're doing GCC or GCHP
get_filename_component(RUNDIR "${RUNDIR}" ABSOLUTE)  # Remove trailing slash if present
get_filename_component(RUNDIR_NAME "${RUNDIR}" NAME)
if("${RUNDIR_NAME}" MATCHES "gchp.*")
    set(IMPL_GUESS "GCHP")
else()
    set(IMPL_GUESS "Classic")
endif()

# Select an implementation
message(STATUS "GEOS-Chem implementation type:")
set_dynamic_option(IMPL "${IMPL_GUESS}"
    LOG IMPL_LOG
    SELECT_EXACTLY 1
    OPTIONS "Classic" "GCHP"
)
dump_log(IMPL_LOG)

# Get the implementations build flags
if("${IMPL}" STREQUAL "GCHP")
    list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/GCHP/CMakeScripts)
    include(GCHPBuildFlags)
else()
    include(GCClassicBuildFlags)
endif()

# Get configuration hash
set(SAFE_CFG_HASHES
    f549694
)

# Sort build definitions (which define build configuration)
set(CFG_SORTED ${GC_DEFINES})
list(SORT CFG_SORTED)

# Get first 7 digits of the hash
string(REPLACE ";" ";" CFG_SORTED "${GC_DEFINES}")
string(SHA1 CFG_HASH "${CFG_SORTED}")
string(SUBSTRING "${CFG_HASH}" 0 7 CFG_HASH)

# Print the configuration hash
message(STATUS "Build configuration hash: ${CFG_HASH}\n")

# Check if the current configuration has been verified
list(FIND SAFE_CFG_HASHES "${CFG_HASH}" SAFE_IDX)
if(${SAFE_IDX} EQUAL -1) 
    message(WARNING "This build configuration, ${CFG_HASH}, has not been validated. Proceed with caution.")
endif()
