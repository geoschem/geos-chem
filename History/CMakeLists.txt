add_static_target(History
    SOURCES
        histcontainer_mod.F90
        histitem_mod.F90
        history_mod.F90
        history_netcdf_mod.F90
        history_util_mod.F90
        metahistcontainer_mod.F90
        metahistitem_mod.F90
        
    PUBLIC_LINK
        KPP
)