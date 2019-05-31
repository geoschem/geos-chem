cmake_minimum_required(VERSION 3.5)
project(GEOS_Chem VERSION 12.3.2 LANGUAGES Fortran C CXX)

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
	gc_message(CONSTANT CMAKE_C_COMPILER DISPLAY_NAME "C")
	gc_message(CONSTANT CMAKE_CXX_COMPILER DISPLAY_NAME "C++")
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

	# Configure for GCClassic
	include(GC-ConfigureClassic)
	configureGCClassic()
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
