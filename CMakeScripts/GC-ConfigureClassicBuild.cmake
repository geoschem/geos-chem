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
find_package(NetCDF REQUIRED)
find_package(OpenMP REQUIRED)

# Set BaseTarget properties
target_compile_options(BaseTarget 
    INTERFACE ${OpenMP_Fortran_FLAGS}
)
target_link_libraries(BaseTarget 
	INTERFACE NetCDF-F ${OpenMP_Fortran_FLAGS}
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

# Inspect SIM and select KPP mech 
set(STANDARD_MECHS
    "standard"
    "benchmark"
    "aciduptake"
    "marinePOA"
    "masscons"
    "TransportTracers"
    "POPs"
    "CH4"
    "tagCH4"
    "tagO3"
    "tagCO"
    "tagHg"
    "CO2"
    "aerosol"
    "Hg"
)
set(TROPCHEM_MECHS
    "tropchem"
    "RRTMG"
    "TOMAS15"
    "TOMAS40"
    "complexSOA"
)
set(SOA_SVPOA_MECHS
    "complexSOA_SVPOA"
)
set(CUSTOM_MECHS
    "custom"
)

inspect_rundir(RUNDIR_SIM 2)

if("${RUNDIR_SIM}" IN_LIST STANDARD_MECHS)
    set(RUNDIR_MECH "Standard")
elseif("${RUNDIR_SIM}" IN_LIST TROPCHEM_MECHS)
    set(RUNDIR_MECH "Tropchem")
elseif("${RUNDIR_SIM}" IN_LIST SOA_SVPOA_MECHS)
    set(RUNDIR_MECH "SOA_SVPOA")
elseif("${RUNDIR_SIM}" IN_LIST CUSTOM_MECHS)
    set(RUNDIR_MECH "custom")
else()
    message(FATAL_ERROR "Unknown simulation type \"${RUNDIR_SIM}\". Cannot determine MECH.")
endif()

# Misc
if("${RUNDIR_SIM}" STREQUAL "masscons")
    set_dynamic_default(GC_DEFINES DEFAULT MASSCONS)
elseif("${RUNDIR_SIM}" MATCHES "TOMAS15")
    set_dynamic_default(GC_DEFINES DEFAULT TOMAS TOMAS15)
    set(NC_DIAG_GUESS "FALSE")
elseif("${RUNDIR_SIM}" MATCHES "TOMAS40")
    set_dynamic_default(GC_DEFINES DEFAULT TOMAS TOMAS40)
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
    set_dynamic_default(GC_DEFINES DEFAULT NESTED_${REGION})
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
    OPTIONS "Standard" "Tropchem" "SOA_SVPOA"
)

# Set reduced grid
set(LAYERS_72_SIMS
    "standard"
    "benchmark"
    "aciduptake"
    "marinePOA"
    "TransportTracers"
    "custom"
)
set(LAYERS_47_SIMS
    "masscons"
    "POPs"
    "CH4"
    "tagCH4"
    "tagO3"
    "tagCO"
    "tagHg"
    "CO2"
    "aerosol"
    "Hg"
    "tropchem"
    "RRTMG"
    "complexSOA"
    "complexSOA_SVPOA"
    "TOMAS15"
    "TOMAS40"
)
if("${RUNDIR_SIM}" IN_LIST LAYERS_72_SIMS)
    set(LAYERS_DEFAULT "72")
elseif("${RUNDIR_SIM}" IN_LIST LAYERS_47_SIMS)
    set(LAYERS_DEFAULT "47")
else()
    message(FATAL_ERROR "Unknown simulation type \"${RUNDIR_SIM}\". Cannot determine LAYERS.")
endif()
set_dynamic_option(LAYERS 
    DEFAULT "${LAYERS_DEFAULT}"
    LOG GENERAL_OPTIONS_LOG
    SELECT_EXACTLY 1
    OPTIONS "47" "72"
)
if("${LAYERS}" STREQUAL "47")
    set_dynamic_default(GC_DEFINES DEFAULT GRIDREDUCED)
endif()

# Build RRTMG?
if("${RUNDIR_SIM}" STREQUAL "RRTMG")
    set(RRTMG_DEFAULT "TRUE")
else()
    set(RRTMG_DEFAULT "FALSE")
endif()
set_dynamic_option(RRTMG 
    DEFAULT ${RRTMG_DEFAULT}
    LOG GENERAL_OPTIONS_LOG
    SELECT_EXACTLY 1
    OPTIONS "TRUE" "FALSE"
)
if(${RRTMG})
    set_dynamic_default(GC_DEFINES DEFAULT "RRTMG")
endif()

# Build GTMM?
set_dynamic_option(GTMM 
    DEFAULT "FALSE"
    LOG GENERAL_OPTIONS_LOG
    SELECT_EXACTLY 1
    OPTIONS "TRUE" "FALSE"
)
if(${GTMM})
    set_dynamic_default(GC_DEFINES DEFAULT "GTMM_Hg")
endif()

if("${RUNDIR_SIM}" STREQUAL "benchmark")
    set(TIMERS_DEFAULT "TRUE")
else()
    set(TIMERS_DEFAULT "FALSE")
endif()
# Build with timers?
set_dynamic_option(TIMERS 
    DEFAULT ${TIMERS_DEFAULT}
    LOG GENERAL_OPTIONS_LOG
    SELECT_EXACTLY 1
    OPTIONS "TRUE" "FALSE"
)
if(${TIMERS})
    set_dynamic_default(GC_DEFINES DEFAULT "USE_TIMERS")
endif()


# Get diagnostics
if("${RUNDIR_SIM}" MATCHES "TOMAS")
    set(DIAG_DEFAULT "BPCH")
else()
    set(DIAG_DEFAULT "NC" "BPCH")
endif()
set_dynamic_option(DIAG
    DEFAULT ${DIAG_DEFAULT}
    OPTIONS "NC" "BPCH"
    LOG GENERAL_OPTIONS_LOG
)

# Add NetCDF diagnostic definitions
if("${DIAG}" MATCHES ".*NC.*")
    set_dynamic_default(GC_DEFINES DEFAULT "NC_DIAG")

    # Read netcdf.inc and search for nf_def_var_deflate
    foreach(NC_INC_DIR ${NETCDF_INCLUDE_DIRS})
        if(EXISTS ${NC_INC_DIR}/netcdf.inc)
            file(READ ${NC_INC_DIR}/netcdf.inc NCINC)
            if("${NCINC}" MATCHES ".*nf_def_var_deflate.*")
                set_dynamic_default(GC_DEFINES DEFAULT "NC_HAS_COMPRESSION")
                break()
            endif()
        endif()
    endforeach()
endif()

# Add BPCH diagnostic definitions
if("${DIAG}" MATCHES ".*BPCH.*")
    set_dynamic_default(GC_DEFINES 
        DEFAULT 
            "BPCH_DIAG" "BPCH_TIMESER" "BPCH_TPBC" 
    )
endif()

# Get flexible precision setting
set_dynamic_option(PREC
    DEFAULT "REAL8"
    SELECT_EXACTLY 1
    OPTIONS "REAL4" "REAL8"
    LOG GENERAL_OPTIONS_LOG
)
if("${PREC}" STREQUAL "REAL8")
    set_dynamic_default(GC_DEFINES DEFAULT "USE_REAL8")
endif()

message(STATUS "General settings:")
dump_log(GENERAL_OPTIONS_LOG)

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
            -cpp                        # Pass through preprocessor before compilation
            -w                          # Disables	all warning messages [attn:Liam this should be removed]
            -auto                       # Causes all local, non-SAVEd variables to be allocated to the run-time stack.
            -noalign                    # Prevents the alignment of data items.
            -convert big_endian         # Specifies the format for unformatted files to be big endian
            -fp-model source            # Prevent any optimizations that would change numerical results
            -mcmodel=medium             # Specifies compilers memory model.
            -shared-intel               # Causes Intel-provided libraries to be linked in dynamically.
            -traceback                  # Generate extra info to provide traceback information when a run time error occurs.
            ${OpenMP_Fortran_FLAGS}
            -DLINUX_IFORT
        LOG RESULTING_DEFINES_LOG
    )
    set(CMAKE_Fortran_FLAGS_RELEASE "-vec-report0 -O2")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g -O0 -check arg_temp_created -debug all -DDEBUG")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    set_dynamic_default(FC_OPTIONS
        DEFAULT
            -cpp 
            -w 
            -std=legacy 
            -fautomatic 
            -fno-align-commons 
            -fconvert=big-endian 
            -fno-range-check
            -mcmodel=medium 
            -fbacktrace 
            -g
            ${OpenMP_Fortran_FLAGS} 
            -DLINUX_GFORTRAN

        LOG RESULTING_DEFINES_LOG
    )
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-loops")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -gdwarf-2 -gstrict-dwarf -O0 -Wall -Wextra -Wconversion -Warray-temporaries -fcheck-array-temporaries")
else()
    message(FATAL_ERROR "${CMAKE_Fortran_COMPILER_ID} Fortran compiler is not currently supported!")
endif()

message(STATUS "Resulting definitions/options:")
dump_log(RESULTING_DEFINES_LOG)

# Set compiler definitions and options in BaseTarget
target_compile_options(BaseTarget 
    INTERFACE 
        ${FC_OPTIONS}
)
unset(FC_OPTIONS)
