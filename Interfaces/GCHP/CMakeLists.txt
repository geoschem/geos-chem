# Interfaces/GCHP/CMakeLists.txt

esma_set_this (OVERRIDE GEOSChem_GridComp)

esma_add_library (${this}
  SRCS Chem_GridCompMod.mapl3.F90
  DEPENDENCIES MAPL GeosCore
  TYPE SHARED)

target_include_directories (${this} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)

mapl_acg (${this} GEOSChem_StateSpecs.acg
  IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS GET_POINTERS DECLARE_POINTERS)

target_compile_definitions (${this} PRIVATE $<$<BOOL:${ADJOINT}>:ADJOINT>)
