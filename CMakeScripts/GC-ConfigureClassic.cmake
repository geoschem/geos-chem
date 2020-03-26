function(configureGCClassic)

    #----------------------------------------------------------------
    # Find OpenMP if we're building a multithreaded executable
    #----------------------------------------------------------------
    gc_pretty_print(SECTION "Threading")
    set(OMP ON CACHE STRING
        "Switch to enable/disable OpenMP threading in GEOS-Chem"
    )
    gc_pretty_print(VARIABLE OMP IS_BOOLEAN)
    if("${OMP}")
       find_package(OpenMP REQUIRED)
       target_compile_options(GEOSChemBuildProperties
           INTERFACE ${OpenMP_Fortran_FLAGS}
       )
       target_link_libraries(GEOSChemBuildProperties
           INTERFACE ${OpenMP_Fortran_FLAGS}
       )
    else()
       target_compile_definitions(GEOSChemBuildProperties
           INTERFACE NO_OMP
       )
    endif()

    #----------------------------------------------------------------
    # Check that GEOS-Chem's version number matches the run directory's version
    #----------------------------------------------------------------
    if(NOT GCCLASSIC_WRAPPER)
        if(EXISTS ${RUNDIR}/Makefile AND NOT "${BUILD_WITHOUT_RUNDIR}")
            # Read ${RUNDIR}/Makefile which has the version number
            file(READ ${RUNDIR}/Makefile RUNDIR_MAKEFILE)
        
            # Pull out the major.minor version
            if(RUNDIR_MAKEFILE MATCHES
        		"VERSION[ \t]*:=[ \t]*([0-9]+\\.[0-9]+)\\.[0-9]+")
                set(RUNDIR_VERSION ${CMAKE_MATCH_1})
            else()
                message(FATAL_ERROR "Failed to determine your run directory's "
                                    "version from ${RUNDIR}/Makefile"
        	    )
            endif()
        
            # Get the major.minor version of GEOS-Chem
            if(PROJECT_VERSION MATCHES "([0-9]+\\.[0-9]+)\\.[0-9]+")
                set(GC_VERSION ${CMAKE_MATCH_1})
            else()
                message(FATAL_ERROR "Internal error. Bad GEOS-Chem version number.")
            endif()
        
            # Throw error if major.minor versions don't match
            if(NOT "${GC_VERSION}" VERSION_EQUAL "${RUNDIR_VERSION}")
                message(FATAL_ERROR
                    "Mismatched version numbers. Your run directory's version "
                    "number is ${RUNDIR_VERSION} but the GEOS-Chem source's "
                    "version number is ${PROJECT_VERSION}"
                )
            endif()
        endif()
    endif()

    #----------------------------------------------------------------
    # Configure the build based on the run directory. 
    #----------------------------------------------------------------
    # Propagate the configuration variables.
    # Define a macro for inspecting the run directory. Inspecting the run
    # directory is how we determine which compiler definitions need to be set.
    macro(inspect_rundir VAR ID)
        if(EXISTS ${RUNDIR}/getRunInfo)
            execute_process(COMMAND perl ${RUNDIR}/getRunInfo ${RUNDIR} ${ID}
                OUTPUT_VARIABLE ${VAR}
                OUTPUT_STRIP_TRAILING_WHITESPACE
            )
        endif()
    endmacro()

    #----------------------------------------------------------------
    # Inspect the run directory to get simulation type
    #----------------------------------------------------------------
    inspect_rundir(RUNDIR_SIM 5)

    #----------------------------------------------------------------
    # Determine the appropriate chemistry mechanism base on the simulation
    #----------------------------------------------------------------
    set(STANDARD_MECHS
        standard      benchmark         aciduptake    marinePOA
        masscons      TransportTracers  POPs          CH4
        tagCH4        tagO3             tagCO
        tagHg         CO2               aerosol
        Hg
        HEMCO # doesn't matter for the HEMCO standalone
    )
    set(TROPCHEM_MECHS
        tropchem      RRTMG     TOMAS15
        TOMAS40       APM       complexSOA
    )
    set(SOA_SVPOA_MECHS
        complexSOA_SVPOA
    )
    set(CUSTOM_MECHS
        custom
    )
    if("${RUNDIR_SIM}" IN_LIST STANDARD_MECHS)
        set(RUNDIR_MECH Standard)
    elseif("${RUNDIR_SIM}" IN_LIST TROPCHEM_MECHS)
        set(RUNDIR_MECH Tropchem)
    elseif("${RUNDIR_SIM}" IN_LIST SOA_SVPOA_MECHS)
        set(RUNDIR_MECH SOA_SVPOA)
    elseif("${RUNDIR_SIM}" IN_LIST CUSTOM_MECHS)
        set(RUNDIR_MECH custom)
    else()
        message(FATAL_ERROR "Unknown simulation type \"${RUNDIR_SIM}\". "
                            "Cannot determine MECH.")
    endif()

    #----------------------------------------------------------------
    # Definitions for specific run directories
    #----------------------------------------------------------------
    set(TOMAS FALSE)
    if("${RUNDIR_SIM}" STREQUAL masscons)
        target_compile_definitions(GEOSChemBuildProperties
	    INTERFACE MASSCONS
	)
    elseif("${RUNDIR_SIM}" MATCHES TOMAS15)
        target_compile_definitions(GEOSChemBuildProperties
	    INTERFACE TOMAS TOMAS15
	)
        set(TOMAS TRUE)
    elseif("${RUNDIR_SIM}" MATCHES TOMAS40)
        target_compile_definitions(GEOSChemBuildProperties
	    INTERFACE TOMAS TOMAS40
	)
        set(TOMAS TRUE)
    endif()

    gc_pretty_print(SECTION "General settings")

    #----------------------------------------------------------------
    # Make MECH an option. This controls which KPP directory is used.
    #----------------------------------------------------------------
    set(MECH "${RUNDIR_MECH}" CACHE STRING "GEOS-Chem's chemistry mechanism")
    # Check that MECH is a valid
    set(VALID_MECHS Standard Tropchem SOA_SVPOA)
    # print MECH to console
    gc_pretty_print(VARIABLE MECH OPTIONS ${VALID_MECHS})
    if(NOT "${MECH}" IN_LIST VALID_MECHS)
        message(FATAL_ERROR "The value of MECH, \"${MECH}\", is an "
			    "invalid chemistry mechanism! Select one "
			    "of: ${VALID_MECHS}."
    )
    endif()

    #----------------------------------------------------------------
    # Turn on bpch diagnostics?
    #----------------------------------------------------------------
    # Define simulations that need bpch diagnostics on
    set(BPCH_ON_SIM RRTMG TOMAS12 TOMAS15 TOMAS30 TOMAS40 Hg tagHg POPs)
    if("${RUNDIR_SIM}" IN_LIST BPCH_ON_SIM)
        set(BPCH_DIAG_DEFAULT ON)
    else()
        set(BPCH_DIAG_DEFAULT OFF)
    endif()
    set(BPCH_DIAG "${BPCH_DIAG_DEFAULT}" CACHE BOOL
    	"Switch to enable GEOS-Chem's bpch diagnostics"
    )
    gc_pretty_print(VARIABLE BPCH_DIAG IS_BOOLEAN)
    if(${BPCH_DIAG})
        target_compile_definitions(GEOSChemBuildProperties
	    INTERFACE BPCH_DIAG
	)
    endif()

    #----------------------------------------------------------------
    # Always set MODEL_CLASSIC when building GEOS-Chem Classic
    #----------------------------------------------------------------
    target_compile_definitions(GEOSChemBuildProperties
	INTERFACE MODEL_CLASSIC
    )

    #----------------------------------------------------------------
    # Always set USE_REAL8. See https://github.com/geoschem/geos-chem/issues/43.
    #----------------------------------------------------------------
    target_compile_definitions(GEOSChemBuildProperties
	INTERFACE USE_REAL8
    )

    #----------------------------------------------------------------
    # Build with timers?
    #----------------------------------------------------------------
    if("${RUNDIR_SIM}" STREQUAL benchmark)
        set(TIMERS_DEFAULT ON)
    else()
        set(TIMERS_DEFAULT OFF)
    endif()
    set(TIMERS "${TIMERS_DEFAULT}" CACHE BOOL
    	"Switch to enable GEOS-Chem's timers"
    )
    gc_pretty_print(VARIABLE TIMERS IS_BOOLEAN)
    # Set USE_TIMERS
    if(${TIMERS})
        target_compile_definitions(GEOSChemBuildProperties
	    INTERFACE USE_TIMERS
        )
    endif()

    gc_pretty_print(SECTION "Components")

    #----------------------------------------------------------------
    # Build APM?
    #----------------------------------------------------------------
    if("${RUNDIR_SIM}" STREQUAL APM)
        set(APM_DEFAULT ON)
    else()
        set(APM_DEFAULT OFF)
    endif()
    set(APM "${APM_DEFAULT}" CACHE BOOL
    	"Switch to build APM as a component of GEOS-Chem"
    )
    gc_pretty_print(VARIABLE APM IS_BOOLEAN)
    if(${APM})
        target_compile_definitions(GEOSChemBuildProperties INTERFACE APM)
    endif()

    #----------------------------------------------------------------
    # Build RRTMG?
    #----------------------------------------------------------------
    if("${RUNDIR_SIM}" STREQUAL RRTMG)
        set(RRTMG_DEFAULT TRUE)
    else()
        set(RRTMG_DEFAULT FALSE)
    endif()
    set(RRTMG "${RRTMG_DEFAULT}" CACHE BOOL
        "Switch to build RRTMG as a component of GEOS-Chem"
    )
    gc_pretty_print(VARIABLE RRTMG IS_BOOLEAN)
    if(${RRTMG})
        target_compile_definitions(GEOSChemBuildProperties INTERFACE RRTMG)
    endif()

    #----------------------------------------------------------------
    # Build GTMM?
    #----------------------------------------------------------------
    set(GTMM OFF CACHE BOOL
        "Switch to build GTMM as a component of GEOS-Chem"
    )
    gc_pretty_print(VARIABLE GTMM IS_BOOLEAN)
    if(${GTMM})
        target_compile_definitions(GEOSChemBuildProperties INTERFACE GTMM_Hg)
    endif()

    #----------------------------------------------------------------
    # Build HEMCO standalone?
    #----------------------------------------------------------------
    if(NOT GCCLASSIC_WRAPPER)
        if("${RUNDIR_SIM}" STREQUAL HEMCO)
            set(HCOSA_DEFAULT TRUE)
        else()
            set(HCOSA_DEFAULT FALSE)
        endif()
        set(HCOSA "${HCOSA_DEFAULT}" CACHE BOOL
            "Switch to build the hemco-standalone (HCOSA) executable"
        )
        gc_pretty_print(VARIABLE HCOSA IS_BOOLEAN)
    endif()

    #----------------------------------------------------------------
    # Build Luo et al wetdep scheme?
    #----------------------------------------------------------------
    set(LUO_WETDEP OFF CACHE BOOL
        "Switch to build the Luo et al (2019) wetdep scheme into GEOS-Chem"
    )
    gc_pretty_print(VARIABLE LUO_WETDEP IS_BOOLEAN)
    if(${LUO_WETDEP})
        target_compile_definitions(GEOSChemBuildProperties
            INTERFACE LUO_WETDEP
        )
    endif()

    #----------------------------------------------------------------
    # Determine which executables should be built
    #----------------------------------------------------------------
    if(NOT GCCLASSIC_WRAPPER)
        set(GCCLASSIC_EXE_TARGETS geos CACHE STRING
            "Executable targets that get built as a part of \"all\""
        )
        if(${HCOSA})
            list(APPEND GCCLASSIC_EXE_TARGETS hemco_standalone)
        endif()
        if(GTMM)
            list(APPEND GCCLASSIC_EXE_TARGETS gtmm)
        endif()
    endif()

    #----------------------------------------------------------------
    # Export the following variables to GEOS-Chem directory's scope
    #----------------------------------------------------------------
    if(NOT GCCLASSIC_WRAPPER)
        set(GCCLASSIC_EXE_TARGETS   ${GCCLASSIC_EXE_TARGETS}    PARENT_SCOPE)
    endif()
    set(GCHP                    FALSE                       PARENT_SCOPE)
    set(MECH                    ${MECH}                     PARENT_SCOPE)
    set(TOMAS                   ${TOMAS}                    PARENT_SCOPE)
    set(APM                     ${APM}                      PARENT_SCOPE)
    set(RRTMG                   ${RRTMG}                    PARENT_SCOPE)
    set(GTMM                    ${GTMM}                     PARENT_SCOPE)
    set(LUO_WETDEP              ${LUO_WETDEP}               PARENT_SCOPE)
    set(RUNDIR                  ${RUNDIR}                   PARENT_SCOPE)
endfunction()
