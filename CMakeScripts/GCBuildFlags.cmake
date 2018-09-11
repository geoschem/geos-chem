#[[ GCBuildFlags.cmake

This module builds the list of compiler options and definitions for 
building GEOS-Chem. Note that compiler options are platform and compiler
dependent, whereas definitions depend on the build configuration (i.e. MET, 
GRID, etc.).

Exported variables:
    GC_DEFINES
    FC_OPTIONS

    RRTMG

]]

# Guard against multiple includes
if(DEFINED GC_BUILD_FLAGS_INCLUDED)
    return()
else()
    set(GC_BUILD_FLAGS_INCLUDED "TRUE")
endif()

# Get run directory
message(STATUS "Run directory setup:")

set_dynamic_default(RUNDIR "<path to run directory>"
    LOG RUNDIR_LOG
)
dump_log(RUNDIR_LOG)

check_path_rules(RUNDIR WARNING_LOG
    EXISTS 
    WRITABLE
    CONTAINS getRunInfo
    IS_OK_RESULT RUNDIR_OK
)

# Error until run directory has been set properly
if(NOT ${RUNDIR_OK})
    stringify_list(WARNING_LOG JOIN "\n" AFTER)
    message(FATAL_ERROR "${WARNING_LOG}")
endif()

# Inspect run directory
include(GCInspectRunDir)


# MET field
set_dynamic_option(MET ${RUNDIR_MET}
    LOG GENERAL_OPTIONS_LOG
    SELECT_EXACTLY 1
    OPTIONS "GEOS_FP" "MERRA2"
)
set_dynamic_default(GC_DEFINES ${MET})

# Check for nested grid
set_dynamic_option(NESTED "${RUNDIR_NESTED}"
    LOG GENERAL_OPTIONS_LOG
    SELECT_EXACTLY 1
    OPTIONS "TRUE" "FALSE"
)
if(${NESTED})
    # Which nested grid?
    set_dynamic_option(REGION "${RUNDIR_REGION}"
        LOG GENERAL_OPTIONS_LOG
        SELECT_EXACTLY 1
        OPTIONS "AS" "CH" "CU" "EU" "NA" 
    )
    set_dynamic_default(GC_DEFINES NESTED NESTED_${REGION})
endif()

# Horizontal grid
if(${NESTED})
    if("${MET}" STREQUAL "MERRA2") # Nested w/ MERRA2 
        set_dynamic_option(GRID "${RUNDIR_GRID}"
            LOG GENERAL_OPTIONS_LOG
            SELECT_EXACTLY 1
            OPTIONS "0.5x0.625"
        )
    else() # Nested w/ GEOS_FP
        set_dynamic_option(GRID "${RUNDIR_GRID}"
            LOG GENERAL_OPTIONS_LOG
            SELECT_EXACTLY 1
            OPTIONS "0.25x0.3125"
        )
    endif()
else() # Not nested
    set_dynamic_option(GRID "${RUNDIR_GRID}"
        LOG GENERAL_OPTIONS_LOG
        SELECT_EXACTLY 1
        OPTIONS "4x5" "2x2.5"
    )
endif()
string(REPLACE "." "" TEMP "GRID${GRID}")
set_dynamic_default(GC_DEFINES ${TEMP})

# Chemistry mechanism
set_dynamic_option(MECH "${RUNDIR_MECH}"
    LOG GENERAL_OPTIONS_LOG
    SELECT_EXACTLY 1
    OPTIONS "Standard" "Tropchem" "SOA_SVPOA"
)

# Build RRTMG?
set_dynamic_option(RRTMG "FALSE"
    LOG GENERAL_OPTIONS_LOG
    SELECT_EXACTLY 1
    OPTIONS "TRUE" "FALSE"
)
if(${RRTMG})
    set_dynamic_default(GC_DEFINES "RRTMG")
endif()

message(STATUS "General settings:")
dump_log(GENERAL_OPTIONS_LOG)


# Get diagnostics
set_dynamic_default(DIAG 
    "BPCH_DIAG" "BPCH_TIMESER" "BPCH_TPBC"

    LOG EXTRA_DEFS_LOG
)
set_dynamic_default(GC_DEFINES ${DIAG})


# Get extra defines
set_dynamic_default(EXTRA 
    "UCX" "USE_REAL8" "USE_TIMERS"
    
    LOG EXTRA_DEFS_LOG
)
set_dynamic_default(GC_DEFINES ${EXTRA})

message(STATUS "Additional definitions:")
dump_log(EXTRA_DEFS_LOG)

# Get resulting GC_DEFINES
string(REPLACE " " ";" GC_DEFINES "${GC_DEFINES}")
set_dynamic_default(GC_DEFINES LOG RESULTING_DEFINES_LOG)

# Get compiler options
set_dynamic_default(FC_OPTIONS
    -fPIC -cpp -w -auto -noalign "-convert big_endian" -O2 -vec-report0 
    "-fp-model source" -openmp -mcmodel=medium -shared-intel -traceback
    -DLINUX_IFORT

    LOG RESULTING_DEFINES_LOG
)

message(STATUS "Resulting definitions/options:")
dump_log(RESULTING_DEFINES_LOG)

# Replace ';' character (delimiting lists) with ' '
string(REPLACE ";" " " FC_OPTIONS "${FC_OPTIONS}")

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