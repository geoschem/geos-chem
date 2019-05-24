cmake_minimum_required(VERSION 3.5)
project(GEOS_Chem VERSION 12.2.0 LANGUAGES Fortran C CXX)

# Set policies
cmake_policy(SET CMP0054 NEW)
cmake_policy(SET CMP0057 NEW)
if(POLICY CMP0074)
	cmake_policy(SET CMP0074 NEW)
endif()

# Add our module path and include our helpers
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/CMakeScripts)
include(GC-Helpers)

# Define the base target (which will store the build configuration properties)
add_library(BaseTarget INTERFACE)

# Set default CMAKE_BUILD_TYPE to Release
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release" 
		CACHE STRING
		"Choose the type of build, options are: Debug Release."
		FORCE
	)
endif()

message(STATUS "Selected compilers:")
message("    C:       ${CMAKE_C_COMPILER}")
message("    C++:     ${CMAKE_CXX_COMPILER}")
message("    Fortran: ${CMAKE_Fortran_COMPILER}")
message(STATUS "Common CMake variables:")
message("  + CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")
set_dynamic_option(CMAKE_BUILD_TYPE
	DEFAULT ${CMAKE_BUILD_TYPE}
	OPTIONS "Release" "Debug"
	LOG CMAKE_VAR_LOG
)
dump_log(CMAKE_VAR_LOG)

# Put all mod files in mod subdirectory
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mod)
target_include_directories(BaseTarget 
	INTERFACE ${PROJECT_BINARY_DIR}/mod
)

# Get the run directory and implementation type 
include(GC-GetRunDir) # sets RUNDIR and IMPL

if("${IMPL}" STREQUAL "Classic")
	include(GC-ConfigureClassicBuild)
elseif("${IMPL}" STREQUAL "GCHP")
    add_subdirectory(GCHP)
endif()

# Add subdirectories
add_subdirectory(KPP)
add_subdirectory(Headers)
add_subdirectory(GeosUtil)
add_subdirectory(NcdfUtil)
add_subdirectory(History)
add_subdirectory(HEMCO)
add_subdirectory(ISOROPIA)
add_subdirectory(GeosRad)
add_subdirectory(GTMM)
add_subdirectory(GeosCore)

