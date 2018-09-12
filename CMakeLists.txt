cmake_minimum_required(VERSION 3.0.2)
project(GEOS_Chem VERSION 12.0.0 LANGUAGES Fortran)

# Add our module path
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeScripts)

# Find NetCDF
find_package(NetCDF REQUIRED COMPONENTS F90)

# Get compiler options and definitions for building GEOS-Chem
include(GCBuildFlags)
set(CMAKE_Fortran_FLAGS ${FC_OPTIONS})

# Add target helpers
include(GCAddTarget)

# Add subdirectories
add_subdirectory(KPP)
add_subdirectory(Headers)
add_subdirectory(GeosUtil)
add_subdirectory(NcdfUtil)

