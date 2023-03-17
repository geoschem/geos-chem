# KPP/Hg/CMakeLists.txt

#-----------------------------------------------------------------------------
# Add libKPP_FirstPass.a
#-----------------------------------------------------------------------------
add_library(KPP_FirstPass STATIC EXCLUDE_FROM_ALL
  gckpp_Precision.F90
  gckpp_Parameters.F90
  gckpp_Monitor.F90
)

# Dependencies
target_link_libraries(KPP_FirstPass
  PUBLIC
  GEOSChemBuildProperties
)

#-----------------------------------------------------------------------------
# Add libKPP.a
#-----------------------------------------------------------------------------
add_library(KPP STATIC EXCLUDE_FROM_ALL
  commonIncludeVars.H
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
  Hg_RateLawFuncs.F90
  Hg_HetStateFuncs.F90
  rateLawUtilFuncs.F90
  #
  # Build these stub modules to avoid compilation errors.
  # These are symbolic links to files in the ../stubs folder.
  #
  stub_aciduptake_DustChemFuncs.F90
  stub_carbon_Funcs.F90
  stub_fullchem_AutoReduceFuncs.F90
  stub_fullchem_HetStateFuncs.F90
  stub_fullchem_SulfurChemFuncs.F90
)

# Dependencies
target_link_libraries(KPP
  PUBLIC
  GeosUtil
)

# Treat REAL as if it were REAL*8
target_compile_options(KPP
  PRIVATE
  ""
  $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:-r8>
  $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-fdefault-real-8 -fdefault-double-8>
)
