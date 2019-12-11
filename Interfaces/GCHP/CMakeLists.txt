# Add the registry
set(MAPL_ACG ${CMAKE_CURRENT_SOURCE_DIR}/../../../../MAPL/MAPL_Base/mapl_acg.pl)
add_subdirectory(Registry)

# Add GCHP targets
add_library(GCHP STATIC
    Chem_GridCompMod.F90
    gigc_chunk_mod.F90
    gigc_historyexports_mod.F90
    gigc_providerservices_mod.F90
)
target_link_libraries(GCHP
	PUBLIC GeosCore
)
add_dependencies(GCHP registry)
target_include_directories(GCHP PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/Registry)
