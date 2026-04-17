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
#if(NOT TARGET GEOSChemBuildProperties)
    add_library(GEOSChemBuildProperties INTERFACE)
#endif()

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
add_subdirectory(GeosRad)
add_subdirectory(GTMM)
add_subdirectory(GeosCore)

#-----------------------------------------------------------------------------
# Build GEOS-Chem Classic driver program if using GC-Classic
#-----------------------------------------------------------------------------
if(MODEL_CLASSIC)
  add_subdirectory(Interfaces/GCClassic)
endif()

#-----------------------------------------------------------------------------
# Build GCHP driver program if using GCHP
#-----------------------------------------------------------------------------
if(MODEL_GCHP)
  add_subdirectory(Interfaces/GCHP)
endif()

#-----------------------------------------------------------------------------
# Additional GEOS-Chem settings if using GCHP (all necessary?)
#-----------------------------------------------------------------------------
if(MODEL_GCHP)

  set(RRTMG OFF CACHE BOOL
    "Switch to build RRTMG as a component of GEOS-Chem"
  )
  set(TOMAS OFF CACHE BOOL
    "Switch to build TOMAS as a component of GEOS-Chem"
  )
  set(LUO_WETDEP OFF CACHE BOOL
    "Switch to build the Luo et al (2020) wetdep scheme into GEOS-Chem"
  )
  set(FASTJX OFF CACHE BOOL
    "Switch to use legacy FAST-JX photolysis in GEOS-Chem"
  )
  set(KPPSA OFF CACHE BOOL
    "Switch to build the KPP-Standalone Box Model"
  )
  
  # Local variables
  set(GTMM                    FALSE)
  set(MECH "fullchem" CACHE STRING "Name of the chemistry mechanism to use")
  
  # Configure build properties for GEOS-Chem
  target_compile_definitions(GEOSChemBuildProperties INTERFACE
    MODEL_GCHPCTM
    MODEL_GCHP
    USE_ESMF
    MAPL_ESMF
    MAPL3
    EXTERNAL_GRID
    NC_HAS_COMPRESSION
    $<$<BOOL:${USE_REAL8}>:USE_REAL8> $<$<BOOL:${RRTMG}>:RRTMG>
    $<$<BOOL:${ADJOINT}>:ADJOINT>
    $<$<BOOL:${REVERSE_OPERATORS}>:REVERSE_OPERATORS> 
    $<$<BOOL:${TOMAS}>:TOMAS>
    $<$<STREQUAL:${TOMAS_BINS},15>:TOMAS15>
    $<$<STREQUAL:${TOMAS_BINS},40>:TOMAS40>
    $<$<BOOL:${LUO_WETDEP}>:LUO_WETDEP>
    $<$<BOOL:${FASTJX}>:FASTJX>
  )
  
  if(NOT MAPL3)
    target_link_libraries(GEOSChemBuildProperties INTERFACE
      $<LINK_ONLY:MAPL.base>
      $<LINK_ONLY:FVdycoreCubed_GridComp>
      $<LINK_ONLY:HCOI_MAPL_ESMF>
      $<$<BOOL:${OMP}>:OpenMP::OpenMP_Fortran>
      MPI::MPI_C MPI::MPI_Fortran
    )
    target_include_directories(GEOSChemBuildProperties INTERFACE
      $<TARGET_PROPERTY:MAPL.base,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:FVdycoreCubed_GridComp,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:HCOI_MAPL_ESMF,INTERFACE_INCLUDE_DIRECTORIES>
      ${CMAKE_Fortran_MODULE_DIRECTORY}
    )
#   GeosCore
  endif()
  if(MAPL3)
    target_link_libraries(GEOSChemBuildProperties INTERFACE
      $<LINK_ONLY:MAPL.base>
      $<LINK_ONLY:FVdycoreCubed_GridComp>
      $<LINK_ONLY:HCOI_MAPL_ESMF>
      $<$<BOOL:${OMP}>:OpenMP::OpenMP_Fortran>
      MPI::MPI_C MPI::MPI_Fortran
      $<LINK_ONLY:mapl3g>
      $<LINK_ONLY:MAPL.generic3g>
      $<LINK_ONLY:MAPL.state>
    )
    target_include_directories(GEOSChemBuildProperties INTERFACE
      $<TARGET_PROPERTY:mapl3g,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:MAPL.generic3g,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:MAPL.state,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:MAPL.base,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:FVdycoreCubed_GridComp,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:HCOI_MAPL_ESMF,INTERFACE_INCLUDE_DIRECTORIES>
      ${CMAKE_Fortran_MODULE_DIRECTORY}
    )
    # Turn off acg for now
    #mapl_acg(GEOSchem_GridComp GEOSchem_StateSpecs.rc
    #  IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS
    #  GET_POINTERS DECLARE_POINTERS
    #  3g
    #)
  endif()

  set(GEOSChem_Fortran_FLAGS_Intel
      -cpp -w -auto -noalign "SHELL:-convert big_endian" "SHELL:-fp-model source"
      -mcmodel=medium -shared-intel -traceback -DLINUX_IFORT
      CACHE STRING "GEOSChem compiler flags for all build types with Intel compilers"
  )
  set(GEOSChem_Fortran_FLAGS_RELEASE_Intel
     -O2
      CACHE STRING "GEOSChem compiler flags for build type release with Intel compilers"
  )
  set(GEOSChem_Fortran_FLAGS_RELWITHDEBINFO_Intel
      -O2
      CACHE STRING "GEOSChem compiler flags for build type relwithdebinfo with Intel compilers"
  )
  set(GEOSChem_Fortran_FLAGS_DEBUG_Intel
      -g -O0 "SHELL:-check arg_temp_created" "SHELL:-debug all" -fpe0 -ftrapuv -check,bounds
      CACHE STRING "GEOSChem compiler flags for build type debug with Intel compilers"
  )
  
  set(GEOSChem_Fortran_FLAGS_GNU
      -cpp -w -std=legacy -fautomatic -fno-align-commons
      -fconvert=big-endian -fno-range-check -mcmodel=medium
      -fbacktrace -g -DLINUX_GFORTRAN -ffree-line-length-none
      CACHE STRING "GEOSChem compiler flags for all build types with GNU compilers"
  )
  set(GEOSChem_Fortran_FLAGS_RELEASE_GNU
     -O3 -funroll-loops
     CACHE STRING "GEOSChem compiler flags for build type release with GNU compilers"
  )
  set(GEOSChem_Fortran_FLAGS_RELWITHDEBINFO_GNU
     -O3 -funroll-loops
     CACHE STRING "GEOSChem compiler flags for build type relwithdebinfo with GNU compilers"
  )
  set(GEOSChem_Fortran_FLAGS_DEBUG_GNU
      -g -gdwarf-2 -gstrict-dwarf -O0 -Wall -Wextra -Wconversion -Warray-temporaries
      -fcheck=array-temps -ffpe-trap=invalid,zero,overflow -finit-real=snan
      -fcheck=bounds -fcheck=pointer
      CACHE STRING "GEOSChem compiler flags for build type debug with GNU compilers"
  )
  set(CMAKE_Fortran_FLAGS "")
  set(CMAKE_Fortran_FLAGS_RELEASE "")
  set(CMAKE_Fortran_FLAGS_Debug "")
  set(CMAKE_Fortran_FLAGS_RelWithDebInfo "")
  
  set(GEOSChem_SUPPORTED_COMPILER_IDS "Intel" "GNU")
  if(NOT CMAKE_Fortran_COMPILER_ID IN_LIST GEOSChem_SUPPORTED_COMPILER_IDS)
     message(FATAL_ERROR "GEOSChem does not support ${CMAKE_Fortran_COMPILER_ID} compilers")
  endif()
  
  #---------------------------------------------------------------------
  # Assign comiler options to build properties
  #---------------------------------------------------------------------
  target_compile_options(GEOSChemBuildProperties
     INTERFACE
     $<$<STREQUAL:${CMAKE_Fortran_COMPILER_ID},Intel>:
        ${GEOSChem_Fortran_FLAGS_Intel}
        $<$<CONFIG:Debug>:${GEOSChem_Fortran_FLAGS_DEBUG_Intel}>
        $<$<CONFIG:RelWithDebInfo>:${GEOSChem_Fortran_FLAGS_RELWITHDEBINFO_Intel}>
        $<$<CONFIG:Release>:${GEOSChem_Fortran_FLAGS_RELEASE_Intel}>
     >
     $<$<STREQUAL:${CMAKE_Fortran_COMPILER_ID},GNU>:
        ${GEOSChem_Fortran_FLAGS_GNU}
        $<$<CONFIG:Debug>:${GEOSChem_Fortran_FLAGS_DEBUG_GNU}>
        $<$<CONFIG:RelWithDebInfo>:${GEOSChem_Fortran_FLAGS_RELWITHDEBINFO_GNU}>
        $<$<CONFIG:Release>:${GEOSChem_Fortran_FLAGS_RELEASE_GNU}>
     >
  )
  
  target_compile_definitions(GEOSChemBuildProperties
      INTERFACE
      $<$<STREQUAL:${CMAKE_Fortran_COMPILER_ID},Intel>:LINUX_IFORT>
      $<$<STREQUAL:${CMAKE_Fortran_COMPILER_ID},GNU>:LINUX_GFORTRAN>
  )
  
  # Print the options that are turned on in GEOS-Chem
  list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_SOURCE_DIR}/CMakeScripts)
  include(GC-Helpers)
  gc_pretty_print(SECTION "Settings")
  gc_pretty_print(VARIABLE MECH OPTIONS "fullchem" "carbon" "custom")
  gc_pretty_print(VARIABLE OMP IS_BOOLEAN)
  gc_pretty_print(VARIABLE USE_REAL8 IS_BOOLEAN)
  gc_pretty_print(VARIABLE APM IS_BOOLEAN)
  gc_pretty_print(VARIABLE TOMAS IS_BOOLEAN)
  gc_pretty_print(VARIABLE RRTMG IS_BOOLEAN)
  gc_pretty_print(VARIABLE GTMM IS_BOOLEAN)
  gc_pretty_print(VARIABLE LUO_WETDEP IS_BOOLEAN)
  gc_pretty_print(VARIABLE FASTJX IS_BOOLEAN)
  gc_pretty_print(VARIABLE KPPSA IS_BOOLEAN)
endif()

