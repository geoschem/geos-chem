# Add libKPPFirstPass.a
add_library(KPPFirstPass_Hg STATIC EXCLUDE_FROM_ALL
      gcHg_Precision.F90
      gcHg_Parameters.F90
      gcHg_Monitor.F90
)
target_link_libraries(KPPFirstPass_Hg
      PUBLIC GEOSChemBuildProperties
)

# Add libKPP.a
add_library(KPP_Hg STATIC EXCLUDE_FROM_ALL
      gcHg_Function.F90
      gcHg_Global.F90
      gcHg_Initialize.F90
      gcHg_Integrator.F90
      gcHg_Jacobian.F90
      gcHg_JacobianSP.F90
      gcHg_LinearAlgebra.F90
      gcHg_Model.F90
      gcHg_Monitor.F90
      gcHg_Parameters.F90
      gcHg_Precision.F90
      gcHg_Rates.F90
      gcHg_Util.F90
)
target_link_libraries(KPP_Hg
      PUBLIC GeosUtil
)

target_compile_options(KPP_Hg PRIVATE
  ""
  # Treat REAL as if it were REAL*8
  $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:-r8>
  $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-fdefault-real-8 -fdefault-double-8>
)
