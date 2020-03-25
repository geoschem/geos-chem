# Add libKPPFirstPass.a
add_library(KPPFirstPass_Standard STATIC EXCLUDE_FROM_ALL
	gckpp_Precision.F90
    gckpp_Parameters.F90
    gckpp_Monitor.F90
)
target_link_libraries(KPPFirstPass_Standard
	PUBLIC GEOSChemBuildProperties
)

# Add libKPP.a
add_library(KPP_Standard STATIC EXCLUDE_FROM_ALL
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
	gckpp_Util.F90
)
target_link_libraries(KPP_Standard
	PUBLIC GeosUtil
)

