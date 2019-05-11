cmake_minimum_required(VERSION 3.0.2)
project(GEOS_Chem VERSION 12.0.0 LANGUAGES Fortran C CXX)

# Set policies
cmake_policy(SET CMP0054 NEW)

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
add_subdirectory(GeosCore)

