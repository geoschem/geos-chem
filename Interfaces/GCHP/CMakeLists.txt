if(MAPL3)

  # MAPL3 (in development): build the GEOS-Chem gridded component as a
  # shared library loaded at run time by MAPL3, with state specs generated
  # from GEOSChem_StateSpecs.acg
  esma_set_this(OVERRIDE GEOSChem_GridComp)

  esma_add_library(${this}
    SRCS Chem_GridCompMod.mapl3.F90 gchp_chunk_mod.F90 gchp_historyexports_mod.F90
    DEPENDENCIES MAPL GeosCore
    TYPE SHARED)

  target_include_directories(${this} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)

  mapl_acg(${this} GEOSChem_StateSpecs.acg
    IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS GET_POINTERS DECLARE_POINTERS)

  target_compile_definitions(${this} PRIVATE $<$<BOOL:${ADJOINT}>:ADJOINT>)

else()

  # MAPL2
  add_subdirectory(Registry)
  add_library(GCHP STATIC
    Chem_GridCompMod.F90
    gchp_chunk_mod.F90
    gchp_historyexports_mod.F90
  )
  target_link_libraries(GCHP
    PUBLIC GeosCore
  )
  add_dependencies(GCHP registry)
  target_include_directories(GCHP PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/Registry)

endif()
