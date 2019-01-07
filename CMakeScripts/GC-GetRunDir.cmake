# Get run directory
message(STATUS "Run directory setup:")
set(RUNDIR_PROMPT "<path to run directory>")
set_dynamic_default(RUNDIR "${RUNDIR_PROMPT}"
    LOG RUNDIR_LOG
	IS_DIRECTORY
)
dump_log(RUNDIR_LOG)
message(STATUS "Bootstrapping ${RUNDIR}")


# Inspect the RUNDIR's name to guess the implementation (GCC or GCHP)
get_filename_component(RUNDIR_NAME "${RUNDIR}" NAME)
if("${RUNDIR_NAME}" MATCHES "gchp.*")
    set(IMPL_GUESS "GCHP")
else()
    set(IMPL_GUESS "Classic")
endif()

# Get the implementation
message(STATUS "GEOS-Chem implementation type:")
set_dynamic_option(IMPL "${IMPL_GUESS}"
    LOG IMPL_LOG
    SELECT_EXACTLY 1
    OPTIONS "Classic" "GCHP"
)
dump_log(IMPL_LOG)
