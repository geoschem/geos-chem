add_library(Hg STATIC EXCLUDE_FROM_ALL
    assignAgeClassToRunningPool.F90 
    assignRanPoolToAgeClass.F90
    CasaRegridModule.F90
    CleanupCASAarrays.F90
    defineArrays.F90
    defineConstants.F90
    doFPARandLAI.F90
    doHerbCarbon.F90
    doHerbCarbonHg.F90
    doHerbivory.F90
    doHgDeposition.F90
    doLatitude.F90
    doLeafRootShedding.F90
    doMaxHg.F90
    doNPP.F90
    doOptimumTemperature.F90
    doPET.F90
    dorestart_mod.F90
    doSoilMoisture.F90
    doTreeCarbon.F90
    doTreeCarbonHg.F90
    getAgeClassBF.F90
    getFireParams.F90
    getFuelWood.F90
    getSoilMoistParams.F90
    getSoilParams.F90
    GTMM_coupled.F90
    HgOutForGEOS.F90
    input_gtmm_mod.F90
    loadCASAinput.F90
    load_GC_data.F90
    loadHgDeposition.F90
    organizeAgeClasses.F90
    processData.F90
    sort_pick_veg.F90
)
target_link_libraries(Hg
    PUBLIC
        Headers
)


if("gtmm" IN_LIST GCCLASSIC_EXE_TARGETS)
    add_executable(gtmm
        GTMM.F90
    )
    target_link_libraries(gtmm
        PUBLIC Hg
    )
    install(TARGETS gtmm
		RUNTIME DESTINATION ${RUNDIR}
	)
endif()
