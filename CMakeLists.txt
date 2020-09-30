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
#-----------------------------------------------------------------------------
# Each subdirectory specifies how it should be built.
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
if(MODEL_GCHPCTM)
    add_subdirectory(Interfaces/GCHP)
endif()
