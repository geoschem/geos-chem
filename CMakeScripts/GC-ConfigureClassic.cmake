function(configureGCClassic)
	# Find NetCDF
	find_package(NetCDF REQUIRED)	
	target_link_libraries(BaseTarget INTERFACE NetCDF-F)

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

    # Configure the compiler options
    if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
        target_compile_options(BaseTarget INTERFACE
            -cpp -w -auto -noalign -convert big_endian -fp-model source -mcmodel=medium
            -shared-intel -traceback -DLINUX_IFORT
        )
        set(CMAKE_Fortran_FLAGS_RELEASE "-O2" PARENT_SCOPE)
        set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -check arg_temp_created -debug all -DDEBUG" PARENT_SCOPE)
    elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
        target_compile_options(BaseTarget INTERFACE
            -cpp -w -std=legacy -fautomatic -fno-align-commons -fconvert=big-endian 
            -fno-range-check -mcmodel=medium -fbacktrace -g -DLINUX_GFORTRAN
        )
        set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-loops" PARENT_SCOPE)
        set(CMAKE_Fortran_FLAGS_DEBUG "-g -gdwarf-2 -gstrict-dwarf -O0 -Wall -Wextra -Wconversion -Warray-temporaries -fcheck-array-temporaries" PARENT_SCOPE)
    else()
        message(FATAL_ERROR "Unknown Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
    endif()
endfunction()
