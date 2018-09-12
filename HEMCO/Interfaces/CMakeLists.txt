add_static_target(HCOI
    SOURCES
        hcoi_esmf_mod.F90
        hcoi_standalone_mod.F90
        hemco_standalone.F90
    
    PUBLIC_LINK
        HCOX
)