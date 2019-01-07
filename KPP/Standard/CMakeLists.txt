# Add libKPPFirstPass.a
add_static_target(KPPFirstPass
    SOURCES 
        gckpp_Precision.F90
        gckpp_Parameters.F90
        gckpp_Monitor.F90
	PUBLIC_LINK
		BaseTarget
)

# Add libKPP.a
add_static_target(KPP
    SOURCES
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
        
    PUBLIC_LINK
        GeosUtil2
)
