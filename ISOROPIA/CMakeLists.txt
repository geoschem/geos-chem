add_static_target(Isoropia
    SOURCES
        isorropiaII_main_mod.F
)
target_compile_options(Isoropia 
    PRIVATE 
        -r8
)
target_compile_definitions(Isoropia
    PUBLIC 
        ${GC_DEFINES}
)

list(APPEND GC_EXTRA_TARGETS "Isoropia")
set(GC_EXTRA_TARGETS "${GC_EXTRA_TARGETS}" PARENT_SCOPE)