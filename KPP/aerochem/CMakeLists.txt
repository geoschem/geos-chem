# Add libKPPFirstPass.a
add_library(KPPFirstPass_aerochem STATIC EXCLUDE_FROM_ALL
      gckpp_Precision.F90
      gckpp_Parameters.F90
      gckpp_Monitor.F90
)
target_link_libraries(KPPFirstPass_aerochem
      PUBLIC GEOSChemBuildProperties
)

# Add libKPP.a
add_library(KPP_aerochem STATIC EXCLUDE_FROM_ALL
      aerochem_RateLawFuncs.F90
      gckpp_AqRates.F90
      gckpp_Function.F90
      gckpp_Global.F90
      gckpp_HetRates.F90
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
      gckpp_Sulfate.F90
      gckpp_Util.F90
      rateLawUtilFuncs.F90
)
target_link_libraries(KPP_aerochem
      PUBLIC GeosUtil
)

target_compile_options(KPP_aerochem PRIVATE
  ""
  # Treat REAL as if it were REAL*8
  $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:-r8>
  $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-fdefault-real-8 -fdefault-double-8>
)
