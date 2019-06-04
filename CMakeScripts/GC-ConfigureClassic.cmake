function(configureGCClassic)
	# Find OpenMP if we're building a multithreaded executable
	gc_message(SECTION "Threading")
	set_dynamic_option(OMP
        DEFAULT "TRUE"
        SELECT_EXACTLY 1
        OPTIONS "TRUE" "FALSE"
        LOG THREADING_LOG
    )
	dump_log(THREADING_LOG)
    if("${OMP}")
		find_package(OpenMP REQUIRED)
		target_compile_options(BaseTarget INTERFACE ${OpenMP_Fortran_FLAGS})
        target_link_libraries(BaseTarget INTERFACE ${OpenMP_Fortran_FLAGS})
    else()
        target_compile_definitions(BaseTarget INTERFACE "NO_OMP")
	endif()

    # Configure the build based on the run directory. Propagate the configuration variables.
    include(GC-ConfigureForClassicRunDirectory)
    configureForClassicRunDirectory()
    set(GCCLASSIC_EXE_TARGETS   ${GCCLASSIC_EXE_TARGETS}    PARENT_SCOPE)
    set(GCHP                    FALSE                       PARENT_SCOPE)
    set(MECH                    ${MECH}                     PARENT_SCOPE)
    set(TOMAS                   ${TOMAS}                    PARENT_SCOPE)
    set(RRTMG                   ${RRTMG}                    PARENT_SCOPE)
    set(GTMM                    ${GTMM}                     PARENT_SCOPE)
    set(RUNDIR                  ${RUNDIR}                   PARENT_SCOPE)
endfunction()
