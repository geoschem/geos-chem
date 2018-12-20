add_static_target(Transport
    SOURCES
        diag_mod.F
        dao_mod.F
        pjc_pfix_mod.F
        tpcore_bc_mod.F
        tpcore_fvdas_mod.F90
        tpcore_window_mod.F90
        pjc_pfix_window_mod.F
        transport_mod.F
        
        # TODO: can this safely be omitted?
        # tendencies_mod.F90
    
    PUBLIC_LINK
        HCOI ${GC_EXTRA_TARGETS}
)
target_compile_options(Transport 
    PRIVATE     
        ""
        $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:-r8>
	$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-fdefault-real-8>
)

add_static_target(GeosCore
    SOURCES
        aero_drydep.F
        aerosol_mod.F
        benchmark_mod.F
        bromocarb_mod.F
        c2h6_mod.F
        carbon_mod.F
        chemistry_mod.F90
        cldice_HBrHOBr_rxn.F
        cleanup.F
        co2_mod.F
        convection_mod.F
        depo_mercury_mod.F
        diag03_mod.F
        diag04_mod.F
        diag1.F
        diag20_mod.F
        diag_2pm.F
        diag3.F
        diag41_mod.F
        diag42_mod.F
        diag48_mod.F
        diag49_mod.F
        diag50_mod.F
        diag51b_mod.F
        diag51_mod.F
        diag53_mod.F
        diag56_mod.F
        diag63_mod.F
        diagnostics_mod.F90
        diag_oh_mod.F
        drydep_mod.F
        dust_mod.F
        emissions_mod.F90
        fast_jx_mod.F
        flexchem_mod.F90
        fvdas_convect_mod.F
        gamap_mod.F
        gc_environment_mod.F90
        geosfp_read_mod.F90
        get_ndep_mod.F
        global_br_mod.F
        global_ch4_mod.F
        gosat_ch4_mod.F
        hcoi_gc_diagn_mod.F90
        hcoi_gc_main_mod.F90
        hco_interface_mod.F90
        initialize.F
        input_mod.F
        isoropiaII_mod.F
        land_mercury_mod.F
        linoz_mod.F
        mapping_mod.F90
        mercury_mod.F
        merra2_read_mod.F90
        mixing_mod.F90
        modis_lai_mod.F90
        ndxx_setup.F
        oasave.F
        ocean_mercury_mod.F
        olson_landmap_mod.F90
        pbl_mix_mod.F
        planeflight_mod.F
        pops_mod.F
        restart_mod.F
        RnPbBe_mod.F
        rpmares_mod.F
        seasalt_mod.F
        set_global_ch4_mod.F90
        set_prof_o3.F
        sfcvmr_mod.F90
        strat_chem_mod.F90
        sulfate_mod.F
        tagged_co_mod.F
        tagged_o3_mod.F
        tccon_ch4_mod.F
        toms_mod.F
        ucx_mod.F
        uvalbedo_mod.F90
        vdiff_mod.F90
        vdiff_pre_mod.F90
        wetscav_mod.F
        YuIMN_Code.F

        # TODO: when to omit these?
        # exchange_mod.F
        # rrtmg_rad_transfer_mod.F
        # tomas_mod.F

    PUBLIC_LINK
        Transport
)

if("${IMPL}" STREQUAL "Classic")
    # Add geos executable target
    add_executable(geos
        main.F
    )
    target_link_libraries(geos
        PUBLIC 
            GeosCore
    )
    install(TARGETS geos
        RUNTIME DESTINATION ${RUNDIR}
    )
endif()

message(STATUS "Initialized executable target geos")
