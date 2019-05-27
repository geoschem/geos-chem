add_library(Isorropia
    isorropiaII_main_mod.F
)
target_link_libraries(Isorropia
	PUBLIC BaseTarget
)
target_compile_options(Isorropia 
    PRIVATE 
        ""
        $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:-r8>
        $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-fdefault-real-8 -fdefault-double-8>
)

list(APPEND GC_EXTRA_TARGETS "Isorropia")
set(GC_EXTRA_TARGETS "${GC_EXTRA_TARGETS}" PARENT_SCOPE)

