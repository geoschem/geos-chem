# Interfaces/GCHP/CMakeLists.txt

# ewl: comment out Registry for mapl3; will implement separately using new mapl_acg
#add_subdirectory(Registry)

add_library (GCHP STATIC
  Chem_GridCompMod.F90
  gchp_chunk_mod.F90
  gchp_historyexports_mod.F90
)
target_link_libraries (GCHP PUBLIC
  GeosCore
)

# ewl: comment out for mapl3; will implement separately using new mapl_acg
#add_dependencies(GCHP registry)
#target_include_directories(GCHP PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/Registry)
