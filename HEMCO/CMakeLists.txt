# HEMCO/CMakeLists.txt

cmake_minimum_required(VERSION 3.5)
project(HEMCO VERSION 3.0.0 LANGUAGES Fortran)

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
# Add CMakeScripts/ to the module path
#---------------------------------------------------------------------
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_SOURCE_DIR}/CMakeScripts)

#---------------------------------------------------------------------
# Declare the HEMCOBuildProperties library
#---------------------------------------------------------------------
# All HEMCO targets depend on HEMCOBuildProperties. This is used to
# control the compiler options and definitions for HEMCO targets
# (via inheritance).
add_library(HEMCOBuildProperties INTERFACE)

#---------------------------------------------------------------------
# Set default compiler options
#---------------------------------------------------------------------
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    target_compile_options(HEMCOBuildProperties
	INTERFACE
        -cpp -w -auto -noalign -convert big_endian -fp-model source -mcmodel=medium
        -shared-intel -traceback -DLINUX_IFORT
    )
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -check arg_temp_created -debug all -DDEBUG")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    target_compile_options(HEMCOBuildProperties
	INTERFACE
        -cpp -w -std=legacy -fautomatic -fno-align-commons -fconvert=big-endian
        -fno-range-check -mcmodel=medium -fbacktrace -g -DLINUX_GFORTRAN
    )
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-loops")
    # same as release for GEOS-Chem
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O3 -funroll-loops")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -gdwarf-2 -gstrict-dwarf -O0 -Wall -Wextra -Wconversion -Warray-temporaries -fcheck-array-temporaries")
else()
    message(FATAL_ERROR "Unknown Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
endif()

#---------------------------------------------------------------------
# Always set USE_REAL8
#---------------------------------------------------------------------
target_compile_definitions(HEMCOBuildProperties
    INTERFACE "USE_REAL8"
)

#---------------------------------------------------------------------
# Set OpenMP setting as cache variable to not override existing definition
#---------------------------------------------------------------------
set(OMP "TRUE" CACHE STRING "Enables/disables HEMCO built with OpenMP")
if("${OMP}")
    find_package(OpenMP REQUIRED)
    target_link_libraries(HEMCOBuildProperties INTERFACE OpenMP::OpenMP_Fortran) # Prefer OpenMP::OpenMP target rather than variables
else()
    target_compile_definitions(HEMCOBuildProperties
        INTERFACE "NO_OMP"
    )
endif()

#---------------------------------------------------------------------
# Put all of HEMCO's mod files in build subdir called mod
#---------------------------------------------------------------------
set(CMAKE_Fortran_MODULE_DIRECTORY ${HEMCO_BINARY_DIR}/mod)
target_include_directories(HEMCOBuildProperties
    INTERFACE ${HEMCO_BINARY_DIR}/mod
)

#---------------------------------------------------------------------
# Append HEMCO's environment variables to CMAKE_PREFIX_PATH
#---------------------------------------------------------------------
list(APPEND CMAKE_PREFIX_PATH
    # Possible NetCDF environment variables
    $ENV{NetCDF_F_ROOT}     $ENV{NetCDF_C_ROOT}     $ENV{NetCDF_ROOT}
    $ENV{NETCDF_F_ROOT}     $ENV{NETCDF_C_ROOT}     $ENV{NETCDF_ROOT}
    $ENV{NetCDF_Fortran_ROOT}
    $ENV{NETCDF_FORTRAN_ROOT}
    
    # Possible HEMCO environmnet variables (still needed? -ewl)
    $ENV{GC_F_BIN} 	    $ENV{GC_BIN}
    $ENV{GC_F_INCLUDE} 	    $ENV{GC_INCLUDE}
    $ENV{GC_F_LIB} 	    $ENV{GC_LIB}
)

#---------------------------------------------------------------------
# Link NetCDF-F to HEMCOBuildProperties
#---------------------------------------------------------------------
find_package(NetCDF REQUIRED)
target_include_directories(HEMCOBuildProperties
    INTERFACE ${NETCDF_INCLUDE_DIRS}
)
target_link_libraries(HEMCOBuildProperties
    INTERFACE ${NETCDF_LIBRARIES}
)

#---------------------------------------------------------------------
# Use the NC_HAS_COMPRESSION def if nf_def_var_deflate is in netcdf.inc
#---------------------------------------------------------------------
if(EXISTS ${NETCDF_F77_INCLUDE_DIR}/netcdf.inc)
    file(READ ${NETCDF_F77_INCLUDE_DIR}/netcdf.inc NCINC)
    if("${NCINC}" MATCHES ".*nf_def_var_deflate.*")
        target_compile_definitions(HEMCOBuildProperties
            INTERFACE "NC_HAS_COMPRESSION"
        )
    endif()
endif()

#---------------------------------------------------------------------
# Print a description of the source code repo's version
# (Needs more work to enable)
#---------------------------------------------------------------------
#get_repo_version(HEMCO_REPO_VERSION ${CMAKE_CURRENT_SOURCE_DIR})
#message(STATUS "HEMCO @ ${HEMCO_REPO_VERSION})

#---------------------------------------------------------------------
# For HEMCO Standalone
#---------------------------------------------------------------------
if(NOT HEMCO_EXTERNAL_CONFIG)

    # This conditional block configures the HEMCO build for HEMCO
    # standalone. It sets HEMCO_EXE_TARGETS and it configures
    # the HEMCOBuildProperties.
    
    # Set CMAKE_BUILD_TYPE to Release by default
    if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE "Release"
            CACHE STRING
            "Set the build type"
            FORCE
    	)
    endif()
    
    # Set the run directory to one directory up
    set(RUNDIR ".." CACHE PATH "Path to your run directory")
    message(STATUS "Bootstrapping  ${RUNDIR}")

    # Make RUNDIR an absolute path
    get_filename_component(RUNDIR "${RUNDIR}"
        ABSOLUTE BASE_DIR "${CMAKE_BINARY_DIR}"
    )
    
    # Configure the run directory, if any.
    set(BUILD_WITHOUT_RUNDIR FALSE)
    if(EXISTS ${RUNDIR}/HEMCO_sa_Spec.rc)
        file(STRINGS ${RUNDIR}/HEMCO_sa_Spec.rc HEMCO_SA)
    elseif(NOT EXISTS RUNDIR)
        set(RUNDIR ${CMAKE_BINARY_DIR})
        set(BUILD_WITHOUT_RUNDIR TRUE)
    else()
        message(FATAL_ERROR "Your run directory doesn't 
                have HEMCO_sa_Spec.rc! Set RUNDIR to a 
                valid run directory.")
    endif()
    
    # Determine which executables should be built
    set(HEMCO_EXE_TARGETS "hemco_standalone"
        CACHE STRING "Executable targets that get built as a part of \"all\""
    )

endif()

#---------------------------------------------------------------------
# For HEMCO as submodule in MAPL/ESMF superproject
#---------------------------------------------------------------------
if(HEMCO_EXTERNAL_CONFIG AND MAPL_ESMF)
    target_link_libraries(HEMCOBuildProperties
        INTERFACE MAPL_Base
    )
    target_compile_definitions(HEMCOBuildProperties
        INTERFACE ESMF_
    )
endif()

#---------------------------------------------------------------------
# Write HEMCOBuildProperties's configuration to a file
#---------------------------------------------------------------------
get_target_property(BT_DEFINITIONS  HEMCOBuildProperties
    INTERFACE_COMPILE_DEFINITIONS
)
get_target_property(BT_OPTIONS      HEMCOBuildProperties
    INTERFACE_COMPILE_OPTIONS
)
get_target_property(BT_LIBRARIES    HEMCOBuildProperties
    INTERFACE_LINK_LIBRARIES
)
get_target_property(BT_INCLUDES     HEMCOBuildProperties
    INTERFACE_INCLUDE_DIRECTORIES
)
file(WRITE ${CMAKE_BINARY_DIR}/HEMCOBuildProperties.txt
    "# This file shows the HEMCOBuildProperties's configuration.\n"
    "\n"
    "HEMCOBuildProperties::INTERFACE_COMPILE_DEFINITIONS:${BT_DEFINITIONS}\n"
    "HEMCOBuildProperties::INTERFACE_COMPILE_OPTIONS:${BT_OPTIONS}\n"
    "HEMCOBuildProperties::INTERFACE_LINK_LIBRARIES:${BT_LIBRARIES}\n"
    "HEMCOBuildProperties::INTERFACE_INCLUDE_DIRECTORIES:${BT_INCLUDES}\n"
    "CMAKE_Fortran_FLAGS_RELEASE:${CMAKE_Fortran_FLAGS_RELEASE}\n"
    "CMAKE_Fortran_FLAGS_DEBUG:${CMAKE_Fortran_FLAGS_DEBUG}\n"
    "CMAKE_Fortran_FLAGS:${CMAKE_Fortran_FLAGS}\n\n"
)

#---------------------------------------------------------------------
# Try to compile a simple program that uses NetCDF-Fortran and OpenMP
#---------------------------------------------------------------------
if(NOT HEMCO_EXTERNAL_CONFIG AND NOT HEMCO_TRY_RUN_PASSED)

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

#---------------------------------------------------------------------
# Add source code directory
#---------------------------------------------------------------------
add_subdirectory(src)
