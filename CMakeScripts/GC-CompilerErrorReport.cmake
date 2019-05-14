
configure_file(
	${CMAKE_SOURCE_DIR}/CMakeScripts/compiler_error_report.sh.in 
	${CMAKE_BINARY_DIR}/compiler_error_report.sh
)
add_custom_command(OUTPUT report.tar.gz
	COMMAND bash compiler_error_report.sh
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
add_custom_target(report 
	DEPENDS report.tar.gz
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
