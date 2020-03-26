# GEOS-Chem repository CMakeLists.txt

cmake_minimum_required(VERSION 3.5)
project(GEOS_Chem VERSION 12.8.0 LANGUAGES Fortran)

#---------------------------------------------------------------------
# Set policies
#---------------------------------------------------------------------
cmake_policy(SET CMP0054 NEW)
cmake_policy(SET CMP0057 NEW)
if(POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif()
if(POLICY CMP0079)
    cmake_policy(SET CMP0079 NEW)
endif()

#---------------------------------------------------------------------
# Add CMakeScripts/ to the module path and import helper functions
#---------------------------------------------------------------------
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_SOURCE_DIR}/CMakeScripts)
include(GC-Helpers)

#---------------------------------------------------------------------
# Declare the GEOSChemBuildProperties library
#---------------------------------------------------------------------
# All GEOS-Chem targets depend on this. This is used to control
# the compiler options and definitions for GEOS-Chem targets
# (via inheritance).
add_library(GEOSChemBuildProperties INTERFACE)

#---------------------------------------------------------------------
# Set default compiler options
#---------------------------------------------------------------------
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    target_compile_options(GEOSChemBuildProperties INTERFACE
        -cpp -w -auto -noalign -convert big_endian -fp-model source -mcmodel=medium
        -shared-intel -traceback -DLINUX_IFORT
    )
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2") # same as release for GEOS-Chem
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -check arg_temp_created -debug all -DDEBUG")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    target_compile_options(GEOSChemBuildProperties INTERFACE
        -cpp -w -std=legacy -fautomatic -fno-align-commons -fconvert=big-endian
        -fno-range-check -mcmodel=medium -fbacktrace -g -DLINUX_GFORTRAN
    )
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-loops")
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O3 -funroll-loops") # same as release for GEOS-Chem
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -gdwarf-2 -gstrict-dwarf -O0 -Wall -Wextra -Wconversion -Warray-temporaries -fcheck-array-temporaries")
else()
    message(FATAL_ERROR "Unknown Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
endif()

#---------------------------------------------------------------------
# Put all of GEOS-Chem's mod files in build subdir called mod
#---------------------------------------------------------------------
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(GEOSChemBuildProperties
    INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/mod
)

#---------------------------------------------------------------------
# Append GEOS-Chem's environment variables to CMAKE_PREFIX_PATH
#---------------------------------------------------------------------
list(APPEND CMAKE_PREFIX_PATH
    # Possible NetCDF environment variables
    $ENV{NetCDF_F_ROOT}     $ENV{NetCDF_C_ROOT} $ENV{NetCDF_ROOT}
    $ENV{NETCDF_F_ROOT}     $ENV{NETCDF_C_ROOT} $ENV{NETCDF_ROOT}
    $ENV{NetCDF_Fortran_ROOT}
    $ENV{NETCDF_FORTRAN_ROOT}
    
    # Possible GEOS-Chem's environmnet variables
    $ENV{GC_F_BIN} 	    $ENV{GC_BIN}
    $ENV{GC_F_INCLUDE} 	    $ENV{GC_INCLUDE}
    $ENV{GC_F_LIB} 	    $ENV{GC_LIB}
)

#---------------------------------------------------------------------
# Link NetCDF-F to GEOSChemBuildProperties
#---------------------------------------------------------------------
find_package(NetCDF REQUIRED)
target_include_directories(GEOSChemBuildProperties INTERFACE
    ${NETCDF_INCLUDE_DIRS}
)
# Not sure if HCOI should be here...
target_link_libraries(GEOSChemBuildProperties INTERFACE
    ${NETCDF_LIBRARIES}
)

#---------------------------------------------------------------------
# Use the NC_HAS_COMPRESSION def if nf_def_var_deflate is in netcdf.inc
#---------------------------------------------------------------------
if(EXISTS ${NETCDF_F77_INCLUDE_DIR}/netcdf.inc)
    file(READ ${NETCDF_F77_INCLUDE_DIR}/netcdf.inc NCINC)
    if("${NCINC}" MATCHES ".*nf_def_var_deflate.*")
        target_compile_definitions(GEOSChemBuildProperties
            INTERFACE "NC_HAS_COMPRESSION"
        )
    endif()
endif()

#---------------------------------------------------------------------
# Print a description of the source code repo's version
#---------------------------------------------------------------------
get_repo_version(GC_REPO_VERSION ${CMAKE_CURRENT_SOURCE_DIR})
message(STATUS "GEOS-Chem @ ${GC_REPO_VERSION}")

#---------------------------------------------------------------------
# For GEOS-Chem Classic only
#---------------------------------------------------------------------
if(NOT GC_EXTERNAL_CONFIG)
    
    # This conditional block configures the GEOS-Chem build
    # for GEOS-Chem Classic. As mentioned above, it sets
    # GCCLASSIC_EXE_TARGETS, RRTMG, GTMM, TOMAS, MECH, and
    # GCHP, and it configures the GEOSChemBuildProperties.
    
    # Set CMAKE_BUILD_TYPE to Release by default
    if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE "Release"
    	    CACHE STRING
            "Set the build type"
    	    FORCE
    	)
    endif()
    
    # Display CMAKE_PREFIX_PATH and CMAKE_BUILD_TYPE
    gc_pretty_print(SECTION "Useful CMake variables")
    gc_pretty_print(VARIABLE CMAKE_PREFIX_PATH)
    gc_pretty_print(VARIABLE CMAKE_BUILD_TYPE)
    
    # Get the run directory
    gc_pretty_print(SECTION "Run directory setup")
    set(RUNDIR ".." CACHE PATH "Path to your run directory")
    gc_pretty_print(VARIABLE RUNDIR)
    message(STATUS "Bootstrapping  ${RUNDIR}")
    # Make RUNDIR an absolute path
    get_filename_component(RUNDIR "${RUNDIR}"
        ABSOLUTE BASE_DIR "${CMAKE_BINARY_DIR}"
    )
    
    # Configure the GEOS-Chem build for GCHP or GC-Classic
    set(BUILD_WITHOUT_RUNDIR FALSE)
    if(EXISTS ${RUNDIR}/input.geos)
        file(STRINGS ${RUNDIR}/input.geos GCHP REGEX ": *gchp_*")
    elseif(EXISTS ${RUNDIR}/getRunInfo)
        # getRunInfo is only in GC-Classic run directories
        set(GCHP FALSE)
    elseif((NOT EXISTS RUNDIR) AND (DEFINED RUNDIR_SIM))
        # special case for building without a run directory
        set(GCHP FALSE)
        set(RUNDIR ${CMAKE_BINARY_DIR})
        set(BUILD_WITHOUT_RUNDIR TRUE)
    else()
        message(FATAL_ERROR "Your run directory doesn't
    		have an input.geos or getRunInfo! Set
    		RUNDIR to a valid run directory."
    	)
    endif()

    # Configure for GCClassic
    include(GC-ConfigureClassic)
    configureGCClassic()
    
endif()

#---------------------------------------------------------------------
# Add all the subdirectories.
#---------------------------------------------------------------------
# Each subdirectory specifies how it should be built.
add_subdirectory(KPP)
add_subdirectory(Headers)
add_subdirectory(GeosUtil)
add_subdirectory(NcdfUtil)
add_subdirectory(History)
add_subdirectory(ObsPack)
add_subdirectory(APM)
if(NOT HEMCO_EXTERNAL_CONFIG)
    add_subdirectory(HEMCO)
endif()
add_subdirectory(ISORROPIA)
add_subdirectory(GeosRad)
add_subdirectory(GTMM)
add_subdirectory(GeosCore)
if(GCHP)
    add_subdirectory(Interfaces/GCHP)
endif()

#---------------------------------------------------------------------
# Write GEOSChemBuildProperties's configuration to a file
#---------------------------------------------------------------------
get_target_property(BT_DEFINITIONS  GEOSChemBuildProperties
    INTERFACE_COMPILE_DEFINITIONS
)
get_target_property(BT_OPTIONS      GEOSChemBuildProperties 
    INTERFACE_COMPILE_OPTIONS
)
get_target_property(BT_LIBRARIES    GEOSChemBuildProperties
    INTERFACE_LINK_LIBRARIES
)
get_target_property(BT_INCLUDES     GEOSChemBuildProperties
    INTERFACE_INCLUDE_DIRECTORIES
)
file(WRITE ${CMAKE_BINARY_DIR}/GEOSChemBuildProperties.txt
    "# This file shows the GEOSChemBuildProperties's configuration.\n"
    "\n"
    "GEOSChemBuildProperties::INTERFACE_COMPILE_DEFINITIONS:${BT_DEFINITIONS}\n"
    "GEOSChemBuildProperties::INTERFACE_COMPILE_OPTIONS:${BT_OPTIONS}\n"
    "GEOSChemBuildProperties::INTERFACE_LINK_LIBRARIES:${BT_LIBRARIES}\n"
    "GEOSChemBuildProperties::INTERFACE_INCLUDE_DIRECTORIES:${BT_INCLUDES}\n"
    "CMAKE_Fortran_FLAGS_RELEASE:${CMAKE_Fortran_FLAGS_RELEASE}\n"
    "CMAKE_Fortran_FLAGS_DEBUG:${CMAKE_Fortran_FLAGS_DEBUG}\n"
    "CMAKE_Fortran_FLAGS:${CMAKE_Fortran_FLAGS}\n\n"
)

#---------------------------------------------------------------------
# Try to compile a simple program that uses NetCDF-Fortran and OpenMP
#---------------------------------------------------------------------
if(NOT GC_EXTERNAL_CONFIG AND NOT GC_TRY_RUN_PASSED)

    # Format definitions with -D prefix
    set(TRY_RUN_DEFINITIONS "")
    foreach(DEF ${BT_DEFINITIONS})
        list(APPEND TRY_RUN_DEFINITIONS "-D${DEF}")
    endforeach()
    
    # Try to compile and run try_compile.F90
    try_run(RUN_FAILED COMPILED_OK
        ${CMAKE_CURRENT_BINARY_DIR}/try_compile                  # binary dir
        ${CMAKE_CURRENT_SOURCE_DIR}/CMakeScripts/try_compile.F90 # test src file
        COMPILE_DEFINITIONS ${TRY_RUN_DEFINITIONS} ${BT_OPTIONS} # compile flags
        LINK_LIBRARIES ${BT_LIBRARIES}	                         # link flags
        CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${BT_INCLUDES}"       # include dirs
    	COMPILE_OUTPUT_VARIABLE COMPILE_OUTPUT
    	RUN_OUTPUT_VARIABLE RUN_OUTPUT
    )

    # Display a warning if its compilation failed
    if(NOT COMPILED_OK)
        if(OMP)
            set(CONDITIONAL_AND_OMP " and OpenMP")
        endif()
        message(WARNING
            "Failed to compile a simple program that uses "
            "NetCDF-Fortran ${CONDITIONAL_AND_OMP}! Could "
            "your NetCDF installation be broken?\nSee "
            "\"FailedCompile.txt\" for more info."
        )
        file(WRITE ${CMAKE_BINARY_DIR}/FailedCompile.txt
            "${COMPILE_OUTPUT}"
        )
    else()
        file(REMOVE ${CMAKE_BINARY_DIR}/FailedCompile.txt)
    endif()
    
    # Display a warning if its execution failed
    if(RUN_FAILED)
        if(OMP)
            set(CONDITIONAL_AND_OMP "and OpenMP ")
        endif()
        message(WARNING
            "A simple program that uses NetCDF-Fortran "
            "${CONDITIONAL_AND_OMP}compiled successfully, "
            "but its execution failed!\n\nSee "
            "\"FailedExecution.txt\" for more info."
        )
        file(WRITE ${CMAKE_BINARY_DIR}/FailedEasyRun.txt
            "${COMPILE_OUTPUT}\n${RUN_OUTPUT}"
        )
    else()
        file(REMOVE ${CMAKE_BINARY_DIR}/FailedEasyRun.txt
            ${CMAKE_BINARY_DIR}/simple_xy.nc
        )
        set(GC_TRY_RUN_PASSED TRUE CACHE INTERNAL
            "try_run passed" FORCE
        )
    endif()
endif()
