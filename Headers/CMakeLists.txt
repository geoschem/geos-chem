add_static_target(Headers
    SOURCES
        charpak_mod.F90
        CMN_DIAG_mod.F
        CMN_FJX_MOD.F
        CMN_O3_mod.F
        CMN_SIZE_mod.F
        diaglist_mod.F90
        errcode_mod.F90
        input_opt_mod.F90
        inquireMod.F90
        physconstants.F
        precision_mod.F
        registry_mod.F90
        registry_params_mod.F90
        species_database_mod.F90
        species_mod.F90
        state_chm_mod.F90
        state_diag_mod.F90
        state_met_mod.F90
    PUBLIC_LINK
        KPPFirstPass
)

if("${IMPL}" STREQUAL "GCHP")
    target_include_directories(Headers
        PUBLIC
	    ${THIRD_PARTY}/include
    )
endif()

target_compile_definitions(Headers
    PUBLIC
        ${GC_DEFINES}
)
