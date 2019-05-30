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

message(STATUS "Selected compilers:")
message("    C:       ${CMAKE_C_COMPILER}")
message("    C++:     ${CMAKE_CXX_COMPILER}")
message("    Fortran: ${CMAKE_Fortran_COMPILER}")
message(STATUS "Common CMake variables:")
message("  + CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")

# Set CMAKE_BUILD_TYPE to Release by default
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release" 
		CACHE STRING
		"Choose the type of build, options are: Debug Release."
		FORCE
	)
endif()

# Make CMAKE_BUILD_TYPE an option. This restricts the possible values
# to Release and Debug.
set_dynamic_option(CMAKE_BUILD_TYPE
	DEFAULT ${CMAKE_BUILD_TYPE}
	OPTIONS "Release" "Debug"
	LOG CMAKE_VAR_LOG
)
dump_log(CMAKE_VAR_LOG)

# Put all mod files in a subdirectory of the build directory called mod
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mod)
target_include_directories(BaseTarget 
	INTERFACE ${PROJECT_BINARY_DIR}/mod
)

# Set the default value of RUNDIR to the parent directory of the build directory
message(STATUS "Run directory setup:")
set(RUNDIR_DEFAULT "..")
set_dynamic_default(RUNDIR DEFAULT "${RUNDIR_DEFAULT}"
    LOG RUNDIR_LOG
	IS_DIRECTORY
)
dump_log(RUNDIR_LOG)
message(STATUS "Bootstrapping ${RUNDIR}")
get_filename_component(RUNDIR "${RUNDIR}" ABSOLUTE)

# Inspect RUNDIR's name to guess the implementation (GCC or GCHP)
get_filename_component(RUNDIR_NAME "${RUNDIR}" NAME)
if("${RUNDIR_NAME}" MATCHES "gchp.*")
    set(IMPL_GUESS "GCHP")
else()
    set(IMPL_GUESS "Classic")
endif()

# Set the default value of the implementation to that guess
message(STATUS "GEOS-Chem implementation type:")
set_dynamic_option(IMPL 
    DEFAULT "${IMPL_GUESS}"
    LOG IMPL_LOG
    SELECT_EXACTLY 1
    OPTIONS "Classic" "GCHP"
)
dump_log(IMPL_LOG)

# Print a description of the source code repo's version
get_repo_version(GC_REPO_VERSION ${CMAKE_SOURCE_DIR})
message(STATUS "GEOS-Chem version: ${GC_REPO_VERSION}")

# The implementation controls who configures BaseTarget
if("${IMPL}" STREQUAL "Classic")
	# Configure BaseTarget for a GC-Classic build
	include(GC-ConfigureClassicBuild)
elseif("${IMPL}" STREQUAL "GCHP")
	# Add GCHP as a subdirectory. The CMake scripts in GCHP setup BaseTarget
    add_subdirectory(GCHP)
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
