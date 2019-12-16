# Build registry
add_custom_command(OUTPUT 
		${CMAKE_CURRENT_BINARY_DIR}/HEMCO_DeclarePointer___.h  
		${CMAKE_CURRENT_BINARY_DIR}/HEMCO_GetPointer___.h  
		${CMAKE_CURRENT_BINARY_DIR}/HEMCO_ImportSpec___.h
		hemco_registry.stamp
	COMMAND ${MAPL_ACG} -v ${CMAKE_CURRENT_SOURCE_DIR}/HEMCO_Registry.rc
	COMMAND cmake -E touch hemco_registry.stamp
	COMMENT "Building HEMCO registry..."
	MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/HEMCO_Registry.rc
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}	
)
add_custom_command(OUTPUT 
		${CMAKE_CURRENT_BINARY_DIR}/GIGCchem_ExportSpec___.h  
		${CMAKE_CURRENT_BINARY_DIR}/GIGCchem_History___.rc    
		${CMAKE_CURRENT_BINARY_DIR}/GIGCchem_InternalSpec___.h
	    ${CMAKE_CURRENT_BINARY_DIR}/GIGCchem_DeclarePointer___.h  
		${CMAKE_CURRENT_BINARY_DIR}/GIGCchem_GetPointer___.h  
		${CMAKE_CURRENT_BINARY_DIR}/GIGCchem_ImportSpec___.h
		chem_registry.stamp
    COMMAND ${MAPL_ACG} -v ${CMAKE_CURRENT_SOURCE_DIR}/Chem_Registry.rc
	COMMAND cmake -E touch chem_registry.stamp
	COMMENT "Building Chem registry..."
	MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/Chem_Registry.rc
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}	
)
add_custom_target(registry DEPENDS hemco_registry.stamp chem_registry.stamp)