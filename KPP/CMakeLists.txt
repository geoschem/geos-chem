# KPP/CMakeLists.txt

#-----------------------------------------------------------------------------
# Build the fullchem mechanism if configured with -DMECH=fullchem
# Also build the KPP-Standalone if configured with -DKPPSA=y
#-----------------------------------------------------------------------------
if("${MECH}" STREQUAL fullchem)
  add_subdirectory(fullchem)
  if(KPPSA)
    add_subdirectory(standalone)
  endif()
endif()

#-----------------------------------------------------------------------------
# Build the carbon mechanism if configured with -DMECH=carbon
#-----------------------------------------------------------------------------
if("${MECH}" STREQUAL carbon)
  add_subdirectory(carbon)
endif()

#-----------------------------------------------------------------------------
# Build the carbon mechanism if configured with -DMECH=carbon_CH4Jacobian
#-----------------------------------------------------------------------------
if("${MECH}" STREQUAL carbon_CH4Jacobian)
  add_subdirectory(carbon_CH4Jacobian)
endif()

#-----------------------------------------------------------------------------
# Build the custom mechanism if configured with -DMECH=custom
# Also build the KPP-Standalone if configured with -DKPPSA=y
#-----------------------------------------------------------------------------
if("${MECH}" STREQUAL custom)
  add_subdirectory(custom)
  if(KPPSA)
    add_subdirectory(standalone)
  endif()
endif()

#-----------------------------------------------------------------------------
# Build the Hg mechanism if configured with -DMECH=Hg
#-----------------------------------------------------------------------------
if("${MECH}" STREQUAL Hg)
  add_subdirectory(Hg)
endif()

# Comment out this option for now (bmy, 9/16/21)
## Build the aciduptake mechanism if configured with -DMECH=aciduptake
#if("${MECH}" STREQUAL aciduptake)
#  add_subdirectory(aciduptake)
#endif()
