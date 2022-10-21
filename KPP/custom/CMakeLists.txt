# KPP/custom/CMakeLists.txt

#----------------------------------------------------------------------------
# Add libKPPFirstPass.a -- custom mechanism
#----------------------------------------------------------------------------
add_library(KPP_FirstPass
  STATIC EXCLUDE_FROM_ALL
  gckpp_Precision.F90
  gckpp_Parameters.F90
  gckpp_Monitor.F90
)

# Define dependencies for libKPP_FirstPass.a
target_link_libraries(KPP_FirstPass PUBLIC GEOSChemBuildProperties)

#----------------------------------------------------------------------------
# Add libKPP.a -- custom mechanism (edit as necessary)
#----------------------------------------------------------------------------
add_library(KPP 
  STATIC EXCLUDE_FROM_ALL
  commonIncludeVars.H
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
  stub_Hg_HetStateFuncs.F90
)

# Define dependencies for libKPP.a
target_link_libraries(KPP PUBLIC GeosUtil)

# When building libKPP.a, treat REAL as if it were REAL*8
target_compile_options(KPP 
  PRIVATE ""
  $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:-r8>
  $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-fdefault-real-8 -fdefault-double-8>
)
