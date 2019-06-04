cmake_minimum_required(VERSION 3.5)
project(GEOS_Chem VERSION 12.4.0 LANGUAGES Fortran)

# Set policies
cmake_policy(SET CMP0054 NEW)
cmake_policy(SET CMP0057 NEW)
if(POLICY CMP0074)
	cmake_policy(SET CMP0074 NEW)
endif()

# Add CMakeScripts/ to the module path and import helper functions
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/CMakeScripts)
include(GC-Helpers)

# Declare the BaseTarget. All GEOS-Chem targets depend on BaseTarget.
# This is used to control the compiler options and definitions for GEOS-Chem
# targets (via inheritance). 
add_library(BaseTarget INTERFACE)

# Put all of GEOS-Chem's mod files in a subdirectory of the build directory called mod
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(BaseTarget 
	INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/mod
)

# Link NetCDF-F to BaseTarget
find_package(NetCDF REQUIRED)
target_link_libraries(BaseTarget INTERFACE NetCDF-F)

# Use the NC_HAS_COMPRESSION definition if nf_def_var_deflate is in netcdf.inc
if(EXISTS ${NETCDF_F77_INCLUDE_DIR}/netcdf.inc)
	file(READ ${NETCDF_F77_INCLUDE_DIR}/netcdf.inc NCINC)
	if("${NCINC}" MATCHES ".*nf_def_var_deflate.*")
		target_compile_definitions(BaseTarget INTERFACE "NC_HAS_COMPRESSION")
	endif()
endif()

# Print a description of the source code repo's version
get_repo_version(GC_REPO_VERSION ${CMAKE_SOURCE_DIR})
message(STATUS "GEOS-Chem version: ${GC_REPO_VERSION}")

if(NOT GC_EXTERNAL_CONFIG)
	# This conditional block configures the GEOS-Chem build for GEOS-Chem 
	# Classic. As mentioned above, it sets GCCLASSIC_EXE_TARGETS, RRTMG, GTMM
	# TOMAS, MECH, and GCHP, and it configures the BaseTarget.

	# Set CMAKE_BUILD_TYPE to Release by default
	if(NOT CMAKE_BUILD_TYPE)
		set(CMAKE_BUILD_TYPE "Release" 
			CACHE STRING
			"Choose the type of build, options are: Debug Release."
			FORCE
		)
	endif()

	# Display the compilers we're using
	gc_message(SECTION "Selected compilers")
	gc_message(CONSTANT CMAKE_Fortran_COMPILER DISPLAY_NAME "Fortran")
	
	# Display CMAKE_PREFIX_PATH and CMAKE_BUILD_TYPE
	gc_message(SECTION "Common CMake variables")
	gc_message(VARIABLE CMAKE_PREFIX_PATH)
	set_dynamic_option(CMAKE_BUILD_TYPE     # making CMAKE_BUILD_TYPE an option restricts it to Release and Debug
		DEFAULT ${CMAKE_BUILD_TYPE}
		OPTIONS "Release" "Debug"
		LOG CMAKE_VAR_LOG
	)
	dump_log(CMAKE_VAR_LOG)

	# Get the run directory
	gc_message(SECTION "Run directory setup")
	set(RUNDIR_DEFAULT "..")
	set_dynamic_default(RUNDIR DEFAULT "${RUNDIR_DEFAULT}"
		LOG RUNDIR_LOG
		IS_DIRECTORY
	)
	dump_log(RUNDIR_LOG)
	message(STATUS "Bootstrapping ${RUNDIR}")
	get_filename_component(RUNDIR "${RUNDIR}" ABSOLUTE) # Make RUNDIR an absolute path

	# Configure the GEOS-Chem build for GCHP or GC-Classic
	if(EXISTS ${RUNDIR}/input.geos)
		file(STRINGS ${RUNDIR}/input.geos GCHP REGEX ": *gchp_*")
		if(GCHP)
			add_subdirectory(GCHP)
		else()
			# Configure for GCClassic
			include(GC-ConfigureClassic)
			configureGCClassic()
		endif()
	else()
		message(FATAL_ERROR "Your run directory doesn't have an input.geos! Set RUNDIR to a valid run directory.")
	endif()

	# Use the default compiler options
	include(GC-DefaultCompilerOptions)
endif()

# Add all the subdirectories. Each subdirectory specifies how it should be built.
add_subdirectory(KPP)
add_subdirectory(Headers)
add_subdirectory(GeosUtil)
add_subdirectory(NcdfUtil)
add_subdirectory(History)
add_subdirectory(ObsPack)
add_subdirectory(HEMCO)
add_subdirectory(ISORROPIA)
add_subdirectory(GeosRad)
add_subdirectory(GTMM)
add_subdirectory(GeosCore)

# Write BaseTarget's configuration to a file
get_target_property(BT_DEFINITIONS  BaseTarget INTERFACE_COMPILE_DEFINITIONS)
get_target_property(BT_OPTIONS      BaseTarget INTERFACE_COMPILE_OPTIONS)
get_target_property(BT_LIBRARIES    BaseTarget INTERFACE_LINK_LIBRARIES)
get_target_property(BT_INCLUDES     BaseTarget INTERFACE_INCLUDE_DIRECTORIES)
file(WRITE ${CMAKE_BINARY_DIR}/BaseTarget.txt
    "# This file shows the BaseTarget's configuration.\n"
    "\n"
    "BaseTarget::INTERFACE_COMPILE_DEFINITIONS:${BT_DEFINITIONS}\n"
    "BaseTarget::INTERFACE_COMPILE_OPTIONS:${BT_OPTIONS}\n"
    "BaseTarget::INTERFACE_LINK_LIBRARIES:${BT_LIBRARIES}\n"
    "BaseTarget::INTERFACE_INCLUDE_DIRECTORIES:${BT_INCLUDES}\n"
    "CMAKE_Fortran_FLAGS_RELEASE:${CMAKE_Fortran_FLAGS_RELEASE}\n"
    "CMAKE_Fortran_FLAGS_DEBUG:${CMAKE_Fortran_FLAGS_DEBUG}\n"
    "CMAKE_Fortran_FLAGS:${CMAKE_Fortran_FLAGS}\n\n"
)