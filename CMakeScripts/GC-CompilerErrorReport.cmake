
set(OG_VALUES "-DCMAKE_COLOR_MAKEFILE=${CMAKE_COLOR_MAKEFILE} -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}")

add_custom_command(OUTPUT configuration.log
    COMMAND ${CMAKE_COMMAND} -DCMAKE_COLOR_MAKEFILE=FALSE -DCMAKE_VERBOSE_MAKEFILE=TRUE . > configuration.log 2>&1
    COMMENT "Logging your configuration..."
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    USES_TERMINAL
)
add_custom_command(OUTPUT build.log
    COMMAND make --no-print-directory VERBOSE=1 > build.log  2>&1
    COMMENT "Logging your build error..."
    DEPENDS configuration.log
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    USES_TERMINAL
)
add_custom_command(OUTPUT returned_variable.stamp
    COMMAND ${CMAKE_COMMAND} ${OG_VALUES} . > returned_variable.stamp 2>&1
    COMMENT "Logging your build error..."
    DEPENDS build.log
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    USES_TERMINAL
)
add_custom_target(report DEPENDS returned_variable.stamp
    COMMAND ${CMAKE_COMMAND} -E tar -czf report.tar.gz build.log configuration.log CMakeCache.txt
    COMMAND ${CMAKE_COMMAND} -E remove build.log configuration.log returned_variable.stamp
    COMMENT "Generating your report..."
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
