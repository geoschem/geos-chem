# Interfaces/GCClassic/CMakeLists.txt

#-----------------------------------------------------------------------------
# Copy the CMake project version to the gc_classic_version.H file
#-----------------------------------------------------------------------------
message("Creating ${CMAKE_CURRENT_SOURCE_DIR}/gc_classic_version.H")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/gc_classic_version.H.in
  ${CMAKE_CURRENT_SOURCE_DIR}/gc_classic_version.H @ONLY)

#-----------------------------------------------------------------------------
# Define libGCClassic.a and specify its dependencies
#-----------------------------------------------------------------------------
add_library(GCClassic
  STATIC EXCLUDE_FROM_ALL
  main.F90
)
target_link_libraries(GCClassic
  PUBLIC
  GeosCore
)
