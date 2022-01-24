# GEOS-Chem repository CMakeLists.txt

cmake_minimum_required(VERSION 3.5)

#-----------------------------------------------------------------------------
# Set policies
#-----------------------------------------------------------------------------
cmake_policy(SET CMP0054 NEW)
cmake_policy(SET CMP0057 NEW)
if(POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif()
if(POLICY CMP0079)
    cmake_policy(SET CMP0079 NEW)
endif()

# Add GEOS-Chem CMake helpers to MODULE_PATH
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_SOURCE_DIR}/CMakeScripts)
include(GC-Helpers)
#-----------------------------------------------------------------------------
# Print header
#-----------------------------------------------------------------------------
get_repo_version(GC_REPO_VERSION ${CMAKE_CURRENT_SOURCE_DIR})
message("=================================================================")
message("GEOS-Chem ${PROJECT_VERSION} (science codebase)")
message("Current status: ${GC_REPO_VERSION}")
message("=================================================================")

#-----------------------------------------------------------------------------
# Declare the GEOSChemBuildProperties library
#
# All GEOS-Chem targets depend on this. This is used to control the compiler
# options and definitions for GEOS-Chem targets (via inheritance).
#
# Optionally, super project can define and configure this target.
#-----------------------------------------------------------------------------
if(NOT TARGET GEOSChemBuildProperties)
    add_library(GEOSChemBuildProperties INTERFACE)
endif()

#-----------------------------------------------------------------------------
# Add all the subdirectories.
# Each subdirectory specifies how it should be built.
#-----------------------------------------------------------------------------
add_subdirectory(KPP)
add_subdirectory(Headers)
add_subdirectory(GeosUtil)
add_subdirectory(NcdfUtil)
add_subdirectory(History)
add_subdirectory(ObsPack)
add_subdirectory(APM)
add_subdirectory(ISORROPIA)
add_subdirectory(GeosRad)
add_subdirectory(GTMM)
add_subdirectory(GeosCore)

#-----------------------------------------------------------------------------
# Build GEOS-Chem Classic driver program (if MODEL_CLASSIC=T)
#-----------------------------------------------------------------------------
if(MODEL_CLASSIC)
  add_subdirectory(Interfaces/GCClassic)
endif()

#-----------------------------------------------------------------------------
# Build GCHPctm driver program (if MODEL_GCHPCTM=T)
#-----------------------------------------------------------------------------
if(MODEL_GCHPCTM)
  add_subdirectory(Interfaces/GCHP)
endif()
