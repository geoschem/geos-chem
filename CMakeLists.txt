cmake_minimum_required(VERSION 3.0.2)
project(GEOS_Chem VERSION 12.0.0 LANGUAGES Fortran CXX)

# Add our module path
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeScripts)
include(GC-Helpers)

# Define the base target which will hold the build configuration
add_library(BaseTarget INTERFACE)

# Put mod files in mod subdirectory
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
