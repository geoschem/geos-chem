# Find OpenMP and make a dependee of BaseTarget
find_package(OpenMP REQUIRED)
target_compile_options(BaseTarget INTERFACE ${OpenMP_Fortran_FLAGS})
target_link_libraries(BaseTarget INTERFACE ${OpenMP_Fortran_FLAGS})

# Macro to call getRunInfo in the run directory
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

if(${MECH} STREQUAL "Tropchem")
    set_dynamic_default(GC_DEFINES GRIDREDUCED)
endif()


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
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    set_dynamic_default(FC_OPTIONS
        -fPIC -cpp -w -auto -noalign -convert big_endian -O2 -vec-report0 
        -fp-model source -mcmodel=medium -shared-intel -traceback -qopenmp
        -DLINUX_IFORT

        LOG RESULTING_DEFINES_LOG
    )
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    set_dynamic_default(FC_OPTIONS
        -cpp -w -std=legacy -fautomatic -fno-align-commons -fconvert=big-endian -fno-range-check -O3 
        -funroll-loops -mcmodel=medium -fbacktrace -g 
        
        -DLINUX_GFORTRAN

        LOG RESULTING_DEFINES_LOG
    )
else()
    message(FATAL_ERROR "${CMAKE_Fortran_COMPILER_ID} Fortran compiler is not currently supported!")
endif()

message(STATUS "Resulting definitions/options:")
dump_log(RESULTING_DEFINES_LOG)


# Set compiler definitions and options in BaseTarget
target_compile_definitions(BaseTarget INTERFACE ${GC_DEFINES})
target_compile_options(BaseTarget INTERFACE ${FC_OPTIONS})
unset(GC_DEFINES)
unset(FC_OPTIONS)
