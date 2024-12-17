# KPP/fullchem/CMakeLists.txt

#----------------------------------------------------------------------------
# Determine KPP integrator name and set corresponding C-preprocessor flags.
# This enables integrator-specific handling in the fullchem mechanism.
#----------------------------------------------------------------------------

# Get integrator name
EXECUTE_PROCESS(
  COMMAND grep \#INTEGRATOR ${CMAKE_CURRENT_SOURCE_DIR}/fullchem.kpp
  OUTPUT_VARIABLE RESULT
)
separate_arguments(SUBSTRINGS UNIX_COMMAND "${RESULT}")
list(GET SUBSTRINGS 1 KPP_INT_NAME)
string(TOLOWER "${KPP_INT_NAME}" KPP_INT_NAME)

# Set C-preprocessor flags
if(KPP_INT_NAME STREQUAL "rosenbrock_autoreduce")
  set(KPP_INT_AUTOREDUCE ON CACHE BOOL "Integrator is rosenbrock_autoreduce?")
  target_compile_definitions(
    GEOSChemBuildProperties INTERFACE KPP_INT_AUTOREDUCE
  )
elseif(KPP_INT_NAME STREQUAL "lsode")
  set(KPP_INT_LSODE ON CACHE BOOL "Integrator is LSODE?")
  target_compile_definitions(
    GEOSChemBuildProperties INTERFACE KPP_INT_LSODE
  )
endif()

# Print results
gc_pretty_print(SECTION "KPP integrator (read from fullchem.kpp)")
gc_pretty_print(VARIABLE KPP_INT_NAME)
gc_pretty_print(VARIABLE KPP_INT_AUTOREDUCE IS_BOOLEAN)
gc_pretty_print(VARIABLE KPP_INT_LSODE IS_BOOLEAN)

#----------------------------------------------------------------------------
# Add libKPPFirstPass.a -- fullchem mechanism
#----------------------------------------------------------------------------
add_library(KPP_FirstPass
  STATIC EXCLUDE_FROM_ALL
  gckpp_Precision.F90
  gckpp_Parameters.F90
  gckpp_Monitor.F90
)

# Dependencies
target_link_libraries(KPP_FirstPass
  PUBLIC
  GEOSChemBuildProperties
)

#----------------------------------------------------------------------------
# Add libKPP.a -- fullchem mechanism
#----------------------------------------------------------------------------
add_library(KPP
  STATIC EXCLUDE_FROM_ALL
  commonIncludeVars.H
  fullchem_AutoReduceFuncs.F90
  fullchem_HetStateFuncs.F90
  fullchem_RateLawFuncs.F90
  fullchem_SulfurChemFuncs.F90
  gckpp_Function.F90
  gckpp_Global.F90
  gckpp_Initialize.F90
  gckpp_Integrator.F90
  gckpp_Jacobian.F90
  gckpp_JacobianSP.F90
  gckpp_LinearAlgebra.F90
  gckpp_Model.F90
  gckpp_Monitor.F90
  gckpp_Parameters.F90
  gckpp_Precision.F90
  gckpp_Rates.F90
  gckpp_Util.F90
  rateLawUtilFuncs.F90
  #
  # Build these stub modules to avoid compilation errors.
  # These are symbolic links to files in the ../stubs folder.
  #
  stub_aciduptake_DustChemFuncs.F90
  stub_carbon_Funcs.F90
  stub_Hg_HetStateFuncs.F90
)

# Dependencies
target_link_libraries(KPP
  PUBLIC
  GeosUtil
)

# Treat REAL as if it were REAL*8
target_compile_options(KPP
  PRIVATE ""
  $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:-r8>
  $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-fdefault-real-8 -fdefault-double-8>
)

# FOR GFORTRAN ONLY: Disable variable tracking for gckpp_Jacobian.F90 in the
# KPP/fullchem mechanism generated with KPP.  This will avoid a compiler
# warning and a restart of the KPP library build.  This only has to be done
# for release types "Release" and "RelWithDebugInfo".
#   -- Bob Yantosca (16 Feb 2023)
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  if(NOT ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Debug"))
    set_source_files_properties(gckpp_Jacobian.F90
      PROPERTIES
      COMPILE_OPTIONS -fno-var-tracking-assignments
    )
  endif()
endif()
