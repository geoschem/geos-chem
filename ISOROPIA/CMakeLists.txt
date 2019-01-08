add_library(Isoropia
    isorropiaII_main_mod.F
)
target_link_libraries(Isoropia
	PUBLIC BaseTarget
)
target_compile_options(Isoropia 
    PRIVATE 
        ""
        $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:-r8>
        $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-fdefault-real-8>
)

list(APPEND GC_EXTRA_TARGETS "Isoropia")
set(GC_EXTRA_TARGETS "${GC_EXTRA_TARGETS}" PARENT_SCOPE)

