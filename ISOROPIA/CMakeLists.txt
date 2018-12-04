add_static_target(Isoropia
    SOURCES
        isorropiaII_main_mod.F
)
target_compile_options(Isoropia 
    PRIVATE 
        ""
        $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:-r8>
        $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-fdefault-real-8>
)
target_compile_definitions(Isoropia
    PUBLIC 
        ${GC_DEFINES}
)

list(APPEND GC_EXTRA_TARGETS "Isoropia")
set(GC_EXTRA_TARGETS "${GC_EXTRA_TARGETS}" PARENT_SCOPE)
