#[[ GC-ConfigureClassicBuild.cmake

This file configures BaseTarget for a GEOS-Chem Classic build. This file does
three things:

    1) Finds dependencies using find_package. For GEOS-Chem Classic these 
       dependencies are:
        a) NetCDF-C and NetCDF-Fortran
        b) OpenMP

    2) Sets the appropriate preprocessor definitions for the run directory.

    3) Sets the default compiler flags.

]]

#[[--------------------------------------------------------------------------]]
#[[     Finding dependencies.                                                ]]
#[[--------------------------------------------------------------------------]]
find_package(NetCDF REQUIRED COMPONENTS F90)
find_package(OpenMP REQUIRED)

# Set BaseTarget properties
target_compile_options(BaseTarget 
    INTERFACE ${OpenMP_Fortran_FLAGS}
)
target_include_directories(BaseTarget
	INTERFACE ${NETCDF_F90_INCLUDE_DIR}
)
target_link_libraries(BaseTarget 
	INTERFACE ${NETCDF_LIBRARIES} ${OpenMP_Fortran_FLAGS}
)

# Print message with the repo's last commit
get_cwd_last_commit_hash(GC_LAST_COMMIT ${CMAKE_SOURCE_DIR})
message(STATUS "GEOS-Chem repository (last commit): ${GC_LAST_COMMIT}")

#[[--------------------------------------------------------------------------]]
#[[     Setting preprocessor definitions.                                    ]]
#[[--------------------------------------------------------------------------]]


#[[     Get defaults for settings by inspecting the run directory.           ]]

# Define a macro to call getRunInfo in the run directory
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


#[[     Settings TUI with defaults from the run directory inspection.        ]]

# MET field
set_dynamic_option(MET 
    DEFAULT ${RUNDIR_MET}
    LOG GENERAL_OPTIONS_LOG
    SELECT_EXACTLY 1
    OPTIONS "GEOS_FP" "MERRA2"
)
set_dynamic_default(GC_DEFINES DEFAULT ${MET})

# Check for nested grid
set_dynamic_option(NESTED 
    DEFAULT "${RUNDIR_NESTED}"
    LOG GENERAL_OPTIONS_LOG
    SELECT_EXACTLY 1
    OPTIONS "TRUE" "FALSE"
)
if(${NESTED})
    # Which nested grid?
    set_dynamic_option(REGION 
        DEFAULT "${RUNDIR_REGION}"
        LOG GENERAL_OPTIONS_LOG
        SELECT_EXACTLY 1
        OPTIONS "AS" "CH" "CU" "EU" "NA" 
    )
    set_dynamic_default(GC_DEFINES DEFAULT NESTED NESTED_${REGION})
endif()

# Horizontal grid
if(${NESTED})
    if("${MET}" STREQUAL "MERRA2") # Nested w/ MERRA2 
        set_dynamic_option(GRID 
            DEFAULT "${RUNDIR_GRID}"
            LOG GENERAL_OPTIONS_LOG
            SELECT_EXACTLY 1
            OPTIONS "0.5x0.625"
        )
    else() # Nested w/ GEOS_FP
        set_dynamic_option(GRID 
            DEFAULT "${RUNDIR_GRID}"
            LOG GENERAL_OPTIONS_LOG
            SELECT_EXACTLY 1
            OPTIONS "0.25x0.3125"
        )
    endif()
else() # Not nested
    set_dynamic_option(GRID 
        DEFAULT "${RUNDIR_GRID}"
        LOG GENERAL_OPTIONS_LOG
        SELECT_EXACTLY 1
        OPTIONS "4x5" "2x2.5"
    )
endif()
string(REPLACE "." "" TEMP "GRID${GRID}")
set_dynamic_default(GC_DEFINES DEFAULT ${TEMP})

# Chemistry mechanism
set_dynamic_option(MECH 
    DEFAULT "${RUNDIR_MECH}"
    LOG GENERAL_OPTIONS_LOG
    SELECT_EXACTLY 1
    OPTIONS "Standard" "Tropchem" "SOA_SVPOA" "benchmark"
)

if(${MECH} STREQUAL "Tropchem")
    set_dynamic_default(GC_DEFINES DEFAULT GRIDREDUCED)
endif()


# Build RRTMG?
set_dynamic_option(RRTMG 
    DEFAULT "FALSE"
    LOG GENERAL_OPTIONS_LOG
    SELECT_EXACTLY 1
    OPTIONS "TRUE" "FALSE"
)
if(${RRTMG})
    set_dynamic_default(GC_DEFINES DEFAULT "RRTMG")
endif()

message(STATUS "General settings:")
dump_log(GENERAL_OPTIONS_LOG)


# Get diagnostics
set_dynamic_default(DIAG
    DEFAULT 
        "BPCH_DIAG" "BPCH_TIMESER" "BPCH_TPBC" 
    LOG EXTRA_DEFS_LOG
)
set_dynamic_default(GC_DEFINES DEFAULT ${DIAG})


# Get extra defines
set_dynamic_default(EXTRA 
    DEFAULT
        "UCX" "USE_REAL8" "USE_TIMERS"
    
    LOG EXTRA_DEFS_LOG
)
set_dynamic_default(GC_DEFINES DEFAULT ${EXTRA})

message(STATUS "Additional definitions:")
dump_log(EXTRA_DEFS_LOG)

# Get resulting GC_DEFINES (overridable)
string(REPLACE " " ";" GC_DEFINES "${GC_DEFINES}")
set_dynamic_default(GC_DEFINES LOG RESULTING_DEFINES_LOG)


#[[     Set resulting defintions on BaseTarget.                              ]]
target_compile_definitions(BaseTarget INTERFACE ${GC_DEFINES})
unset(GC_DEFINES)


#[[--------------------------------------------------------------------------]]
#[[     Setting default compiler options.                                    ]]
#[[--------------------------------------------------------------------------]]

if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    set_dynamic_default(FC_OPTIONS
        DEFAULT
            -fPIC -cpp -w -auto -noalign -convert big_endian -vec-report0 
            -fp-model source -mcmodel=medium -shared-intel -traceback -qopenmp
            -DLINUX_IFORT

        LOG RESULTING_DEFINES_LOG
    )
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    set_dynamic_default(FC_OPTIONS
        DEFAULT
            -cpp -w -std=legacy -fautomatic -fno-align-commons -fconvert=big-endian -fno-range-check 
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
target_compile_options(BaseTarget 
    INTERFACE 
        ${FC_OPTIONS}

        # Debug flags
        $<$<AND:$<CXX_COMPILER_ID:Intel>,$<CONFIG:DEBUG>>:-g -O0>
        $<$<AND:$<CXX_COMPILER_ID:GNU>,  $<CONFIG:DEBUG>>:-g -Og>

        # Release flags
        $<$<AND:$<CXX_COMPILER_ID:Intel>,$<CONFIG:DEBUG>>:-O2>
        $<$<AND:$<CXX_COMPILER_ID:GNU>,  $<CONFIG:DEBUG>>:-O3>
)
unset(FC_OPTIONS)
