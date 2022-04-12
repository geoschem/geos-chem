# KPP/CCycle/CMakeLists.txt

#----------------------------------------------------------------------------
# Add libKPPFirstPass.a -- fullchem mechanism
#----------------------------------------------------------------------------
add_library(CCYCLE_FirstPass
  STATIC EXCLUDE_FROM_ALL
  ccycle_Precision.F90
  ccycle_Parameters.F90
  ccycle_Monitor.F90
)

# Define dependencies for libKPP_FirstPass.a
target_link_libraries(CCYCLE_FirstPass PUBLIC GEOSChemBuildProperties)

#----------------------------------------------------------------------------
# Add libKPP.a -- fullchem mechanism
#----------------------------------------------------------------------------
add_library(CCYCLE 
  STATIC EXCLUDE_FROM_ALL
  ccycle_Function.F90
  ccycle_Global.F90
  ccycle_Initialize.F90
  ccycle_Integrator.F90
  ccycle_Jacobian.F90
  ccycle_JacobianSP.F90
  ccycle_LinearAlgebra.F90
  ccycle_Model.F90
  ccycle_Monitor.F90
  ccycle_Parameters.F90
  ccycle_Precision.F90
  ccycle_Rates.F90
  ccycle_Util.F90
)

# Add dependencies
target_link_libraries(CCYCLE PUBLIC GeosUtil)

# When building libKPP.a, treat REAL as if it were REAL*8
target_compile_options(CCYCLE 
  PRIVATE ""
  $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:-r8>
  $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-fdefault-real-8 -fdefault-double-8>
)
