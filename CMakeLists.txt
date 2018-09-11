cmake_minimum_required(VERSION 3.0.2)
project(GEOS_Chem VERSION 12.0.0 LANGUAGES Fortran)

# Add our module path
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeScripts)

# Include logging and setting helpers
include(GCHelpers)

# Get compiler options and definitions for building GEOS-Chem
include(GCBuildFlags)
include(GCBuildFlags)
