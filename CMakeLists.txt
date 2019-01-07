cmake_minimum_required(VERSION 3.0.2)
project(GEOS_Chem VERSION 12.0.0 LANGUAGES Fortran)

# Add our module path
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeScripts)
include(GC-Helpers)

# Find NetCDF
find_package(NetCDF REQUIRED COMPONENTS F90)

# Define the base target which will hold the build configuration
add_library(BaseTarget INTERFACE)

# Get the run directory and implementation type
include(GC-GetRunDir) # sets RUNDIR and IMPL

if("${IMPL}" STREQUAL "Classic")
	include(GC-ConfigureClassicBuild)
endif()

# Add target helpers
include(GCAddTarget)

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
if("${IMPL}" STREQUAL "GCHP")
    add_subdirectory(GCHP)
endif()

