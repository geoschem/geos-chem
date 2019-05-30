# Find NetCDF on the local machine. Make NetCDF-F a dependency of BaseTarget.
find_package(NetCDF REQUIRED)
target_link_libraries(BaseTarget 
	INTERFACE NetCDF-F
)

# Define a macro for inspecting the run directory. Inspecting the run
# directory is how we determine which compiler definitions need to be set.
macro(inspect_rundir VAR ID)
    execute_process(COMMAND perl ${RUNDIR}/getRunInfo ${RUNDIR} ${ID}
        OUTPUT_VARIABLE ${VAR}
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
endmacro()

# Inspect the run directory to get the met field type and grid resolution
inspect_rundir(RUNDIR_MET 0)
if("${RUNDIR_MET}" STREQUAL "geosfp")
    set(RUNDIR_MET "GEOS_FP")
elseif("${RUNDIR_MET}" STREQUAL "merra2")
    set(RUNDIR_MET "MERRA2")
endif()
inspect_rundir(RUNDIR_GRID 1)
if("${RUNDIR_GRID}" STREQUAL "2x25")
    set(RUNDIR_GRID "2x2.5")
elseif("${RUNDIR_GRID}" STREQUAL "05x0625")
    set(RUNDIR_GRID "0.5x0.625")
elseif("${RUNDIR_GRID}" STREQUAL "025x03125")
    set(RUNDIR_GRID "0.25x0.3125")
endif()

# Inspect the run directory to get simulation type
inspect_rundir(RUNDIR_SIM 2)

# Determine the appropriate chemistry mechanism base on the simulation
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
    "HEMCO" # doesn't matter for the HEMCO standalone
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

# Definitions for specific run directories
if("${RUNDIR_SIM}" STREQUAL "masscons")
    set_dynamic_default(GC_DEFINES DEFAULT MASSCONS)
elseif("${RUNDIR_SIM}" MATCHES "TOMAS15")
    set_dynamic_default(GC_DEFINES DEFAULT TOMAS TOMAS15)
    set(NC_DIAG_GUESS "FALSE")
elseif("${RUNDIR_SIM}" MATCHES "TOMAS40")
    set_dynamic_default(GC_DEFINES DEFAULT TOMAS TOMAS40)
endif()

# Inspect the run directory to determine if it's a nested simulation
inspect_rundir(RUNDIR_REGION 3)
if("${RUNDIR_REGION}" STREQUAL "n")
    set(RUNDIR_NESTED "FALSE")
    unset(RUNDIR_REGION)
else()
    set(RUNDIR_NESTED "TRUE")
    string(TOUPPER "${RUNDIR_REGION}" RUNDIR_REGION)
endif()

# Make MET an option and set the appropriate definition
set_dynamic_option(MET 
    DEFAULT ${RUNDIR_MET}
    LOG GENERAL_OPTIONS_LOG
    SELECT_EXACTLY 1
    OPTIONS "GEOS_FP" "MERRA2"
)
set_dynamic_default(GC_DEFINES DEFAULT ${MET})

# Make NESTED an option and set the appropriate definitions
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

# Make GRID an option with different options based on MET and NESTED, and
# set the appropriate definition
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

# Make MECH an option. This controls which KPP directory is used.
set_dynamic_option(MECH 
    DEFAULT "${RUNDIR_MECH}"
    LOG GENERAL_OPTIONS_LOG
    SELECT_EXACTLY 1
    OPTIONS "Standard" "Tropchem" "SOA_SVPOA"
)

# Make LAYERS an option and set the appropriate definitions. Determine the 
# default value based on RUNDIR_SIM.
set(LAYERS_72_SIMS
    "standard"
    "benchmark"
    "aciduptake"
    "marinePOA"
    "TransportTracers"
    "custom"
    "HEMCO" # doesn't matter for the HEMCO standalone
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

# Build hemco_standalone?
if("${RUNDIR_SIM}" STREQUAL "HEMCO")
    set(HCOSA_DEFAULT "TRUE")
else()
    set(HCOSA_DEFAULT "FALSE")
endif()
set_dynamic_option(HCOSA 
    DEFAULT "${HCOSA_DEFAULT}"
    LOG GENERAL_OPTIONS_LOG
    SELECT_EXACTLY 1
    OPTIONS "TRUE" "FALSE"
)

# Build with timers?
if("${RUNDIR_SIM}" STREQUAL "benchmark")
    set(TIMERS_DEFAULT "TRUE")
else()
    set(TIMERS_DEFAULT "FALSE")
endif()
set_dynamic_option(TIMERS 
    DEFAULT ${TIMERS_DEFAULT}
    LOG GENERAL_OPTIONS_LOG
    SELECT_EXACTLY 1
    OPTIONS "TRUE" "FALSE"
)
if(${TIMERS})
    set_dynamic_default(GC_DEFINES DEFAULT "USE_TIMERS")
endif()


# Build with BPCH diagnostics?
set_dynamic_option(BPCH_DIAG 
    DEFAULT "TRUE"
    OPTIONS "TRUE" "FALSE"
    LOG GENERAL_OPTIONS_LOG
)
if(${BPCH_DIAG})
    set_dynamic_default(GC_DEFINES DEFAULT "BPCH_DIAG" "BPCH_TIMESER" "BPCH_TPBC")
endif()

# Use the NC_HAS_COMPRESSION definition if nf_def_var_deflate is in netcdf.inc
if(EXISTS ${NETCDF_F77_INCLUDE_DIR}/netcdf.inc)
    file(READ ${NETCDF_F77_INCLUDE_DIR}/netcdf.inc NCINC)
    if("${NCINC}" MATCHES ".*nf_def_var_deflate.*")
        set_dynamic_default(GC_DEFINES DEFAULT "NC_HAS_COMPRESSION")
    endif()
endif()

# Make an option for controlling the flexible precision. Set the appropriate
# definition
set_dynamic_option(PREC
    DEFAULT "REAL8"
    SELECT_EXACTLY 1
    OPTIONS "REAL4" "REAL8"
    LOG GENERAL_OPTIONS_LOG
)
if("${PREC}" STREQUAL "REAL8")
    set_dynamic_default(GC_DEFINES DEFAULT "USE_REAL8")
endif()

# Build a single-threaded or multi-threaded executable?
set_dynamic_option(OMP
    DEFAULT "TRUE"
    SELECT_EXACTLY 1
    OPTIONS "TRUE" "FALSE"
    LOG GENERAL_OPTIONS_LOG
)
if("${OMP}")
    find_package(OpenMP REQUIRED)
    set(OMP_FORTRAN_FLAGS "${OpenMP_Fortran_FLAGS}")
    target_link_libraries(BaseTarget INTERFACE ${OpenMP_Fortran_FLAGS})
else()
    set(OMP_FORTRAN_FLAGS "")
    set_dynamic_default(GC_DEFINES DEFAULT "NO_OMP")
endif()

message(STATUS "General settings:")
dump_log(GENERAL_OPTIONS_LOG)

# By using set_dynamic_default, the resulting GC_DEFINES is overwritable
string(REPLACE " " ";" GC_DEFINES "${GC_DEFINES}")
set_dynamic_default(GC_DEFINES LOG RESULTING_DEFINES_LOG)

# Set the definitions for the BaseTarget
target_compile_definitions(BaseTarget INTERFACE ${GC_DEFINES})
unset(GC_DEFINES)


# Set the Fortran compiler options based on the compiler's family
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
            ${OMP_FORTRAN_FLAGS}
            -DLINUX_IFORT
        LOG RESULTING_DEFINES_LOG
    )
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
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
            ${OMP_FORTRAN_FLAGS}
            -DLINUX_GFORTRAN

        LOG RESULTING_DEFINES_LOG
    )
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-loops")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -gdwarf-2 -gstrict-dwarf -O0 -Wall -Wextra -Wconversion -Warray-temporaries -fcheck-array-temporaries")
else()
    message(FATAL_ERROR "${CMAKE_Fortran_COMPILER_ID} Fortran compiler is currently not supported!")
endif()

message(STATUS "Resulting definitions/options:")
dump_log(RESULTING_DEFINES_LOG)

# Set compiler options for the BaseTarget
target_compile_options(BaseTarget INTERFACE ${FC_OPTIONS})
unset(FC_OPTIONS)