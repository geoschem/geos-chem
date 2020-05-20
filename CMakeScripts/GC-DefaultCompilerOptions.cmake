if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    target_compile_options(GEOSChemBuildProperties
	INTERFACE
        -cpp -w -auto -noalign -convert big_endian -fp-model source -mcmodel=medium
        -shared-intel -traceback -DLINUX_IFORT
    )
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2") # same as release for GEOS-Chem
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -check arg_temp_created -debug all -fpe0 -ftrapuv -check bounds -traceback")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    target_compile_options(GEOSChemBuildProperties
	INTERFACE
        -cpp -w -std=legacy -fautomatic -fno-align-commons -fconvert=big-endian
        -fno-range-check -mcmodel=medium -fbacktrace -g -DLINUX_GFORTRAN
    )
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-loops")
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O3 -funroll-loops") # same as release for GEOS-Chem
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -gdwarf-2 -gstrict-dwarf -O0 -Wall -Wextra -Wconversion -Warray-temporaries -fcheck-array-temporaries -ffpe-trap=invalid,zero,overflow -finit-real=snan -fbounds-check -fbacktrace")
else()
    message(FATAL_ERROR "Unknown Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
endif()
