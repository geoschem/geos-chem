
#[[ stringify_list

Stringify a list of strings.

Usage:
    stringify_list(<list> 
        [PRINT] 
        [LINE_LENGTH <length>] 
        [HIGHLIGHT <keyword> ...] 
        [JOIN <token> ...] 
        [AFTER]
    )

Options:
    PRINT           Print the stringified list to console. Highlighted text 
                    will also be colorized.
    
    LINE_LENGTH     When JOINing a list of string, the resulting lines will
                    be limited to <length> characters. The resulting <list>
                    will be a list of lines that can then be JOINed with 
                    newlines.
    
    HIGHLIGHT       A list of keywords to highlight.

    JOIN            A list of tokens that will be used sequentially to join
                    list items. The last token will be used to join all 
                    remaining items.

    AFTER           Place the JOIN tokens after the item, rather than before.

]]
function(stringify_list LIST)
    cmake_parse_arguments(BETTER
        "PRINT;AFTER" 
        "LINE_LENGTH" 
        "HIGHLIGHT;JOIN" 
        ${ARGN}
    )

    if(NOT DEFINED BETTER_LINE_LENGTH)
        set(BETTER_LINE_LENGTH 1000) # Arbitrary big number
    endif()

    set(STR ${${LIST}})
    
    # Limit joined line length
    if(DEFINED BETTER_JOIN)
        set(TEMP "")
        set(CUR_LEN "0")

        set(JOIN_IDX "0 - 1")
        list(LENGTH BETTER_JOIN JOIN_LEN)
        
        foreach(ITEM ${STR})
            # Get the join token    
            math(EXPR JOIN_IDX "${JOIN_IDX} + 1")
            if(${JOIN_IDX} LESS ${JOIN_LEN})
                list(GET BETTER_JOIN "${JOIN_IDX}" JOIN_TOKEN)
            endif()
            string(LENGTH "${JOIN_TOKEN}" SEP_LEN)

            # If a line length was 
            string(LENGTH "${ITEM}" WORD_LEN)
            math(EXPR POST_LEN "${WORD_LEN} + ${CUR_LEN} + ${SEP_LEN}")
            if("${POST_LEN}" LESS "${BETTER_LINE_LENGTH}")
                if(${BETTER_AFTER})
                    set(TEMP "${TEMP}${ITEM}${JOIN_TOKEN}")
                else()
                    set(TEMP "${TEMP}${JOIN_TOKEN}${ITEM}")
                endif()
                set(CUR_LEN "${POST_LEN}")
            else()
                if(${BETTER_AFTER})
                    set(TEMP "${TEMP};${ITEM}${JOIN_TOKEN}")
                else()
                    set(TEMP "${TEMP};${JOIN_TOKEN}${ITEM}")
                endif()
                set(CUR_LEN "0")
                math(EXPR CUR_LEN "${SEP_LEN} + ${WORD_LEN}")
            endif()
        endforeach()

        set(STR "${TEMP}")
    endif()

    # Highlight selected words
    if(DEFINED BETTER_HIGHLIGHT)
        foreach(KEYWORD ${BETTER_HIGHLIGHT})
            string(REPLACE "${KEYWORD}" "[${KEYWORD}]" STR "${STR}")
        endforeach()
    endif()
    
    if(${BETTER_PRINT})
        string(ASCII 27 Esc)
        if(${CMAKE_COLOR_MAKEFILE})
            string(REGEX REPLACE "\\[([a-zA-Z0-9_\\.]+)\\]" "${Esc}[32m\\1${Esc}[m" COLORIZED "${STR}")
        else()
            set(COLORIZED "${STR}")
        endif()
	string(REGEX REPLACE "\n$" "" COLORIZED "${COLORIZED}")
        message("${COLORIZED}")
    endif()

    
    # Export the new string
    set(${LIST} "${STR}" PARENT_SCOPE)
endfunction()

#[[ dump_log

Prints a list of logged lines to the console.

Usage:
    dump_log(LOG)

]]
macro(dump_log LOG)
    if(DEFINED ${LOG})
        stringify_list(${LOG} 
            JOIN "\n" AFTER
            PRINT
        )
    endif()
endmacro()

#[[ set_dynamic_default

Sets/modifies VAR iff the user has not overwritten it. This function extends
the basic functionality of CMake's set, such that modifications to VAR will 
only be made iff the variable isn't cached.

Usage:
    set_dynamic_default(VAR <default value> ... 
        [LOG <log>]
    )

Options:
    LOG     The value of VAR will be logged to <log>.

]]
function(set_dynamic_default VAR)
    cmake_parse_arguments(SDD
		"IS_DIRECTORY"
        "LOG"
        "DEFAULT"
        ${ARGN}
    )

    # Set flag indicating if ${VAR} can be set
    if(NOT DEFINED ${VAR} AND NOT DEFINED ${VAR}_IS_MUTABLE)
        set(${VAR}_IS_MUTABLE "true")
        set(${VAR}_IS_MUTABLE "true" PARENT_SCOPE)
    endif()

    # If ${VAR} is mutable, export it
    if(${${VAR}_IS_MUTABLE})
        set(${VAR} ${${VAR}} ${SDD_DEFAULT})
        set(${VAR} ${${VAR}} PARENT_SCOPE)
    endif()

    # Log if requested
    if(DEFINED SDD_LOG)
        set(STR "${${VAR}}")

        # Split list with "  
        stringify_list(STR 
            JOIN "  " 
            LINE_LENGTH 60
        )
        # Wrap lines
        stringify_list(STR 
            JOIN "  + ${VAR}:\t" "\n  ...       \t"
        )
        list(APPEND ${SDD_LOG} "${STR}")
        set(${SDD_LOG} ${${SDD_LOG}} PARENT_SCOPE)
    endif()

	# if IS_DIRECTORY, check that the specified path is valid.
	if(SDD_IS_DIRECTORY)
        # If not absolute, force it to
        if(NOT IS_ABSOLUTE "${${VAR}}")
            set(${VAR} ${CMAKE_BINARY_DIR}/${${VAR}})
            set(${VAR} ${${VAR}} PARENT_SCOPE)
		endif()

		if(NOT EXISTS "${${VAR}}")
			dump_log(${SDD_LOG})
			message(FATAL_ERROR "Invalid ${VAR}. ${${VAR}} does not exist!")
		endif()
	endif()
endfunction()

#[[ set_dynamic_option

Sets variables whose values must be selected from a set of possible values.
Setting and modifying behaviour is the same as set_dynamic_default.

Usage:
    set_dynamic_option(VAR <default value> ... 
        [LOG <log>]
        [SELECT_AT_LEAST <num>]
        [SELECT_AT_MOST <num>]
        [SELECT_EXACTLY <num>]
        [OPTIONS <option> ...]
    )

Options:
    LOG             The value of VAR will be logged to <log>.

    SELECT_AT_MOST  Assert that at least <num> options are selected.

    SELECT_AT_MOST  Assert that at most <num> options are selected.

    SELECT_EXACTLY  Assert that exactly <num> options are selected.

    OPTIONS         The list of valid options.
    
]]
function(set_dynamic_option VAR)
    cmake_parse_arguments(SDO  
        ""
        "SELECT_AT_LEAST;SELECT_AT_MOST;SELECT_EXACTLY;LOG" 
        "OPTIONS;DEFAULT" 
        ${ARGN}
    )
    
    # Set/get value
    set_dynamic_default(${VAR} DEFAULT ${SDO_DEFAULT})
    if(DEFINED SDO_LOG)
        set(STR "${SDO_OPTIONS}")

        # Split list with "  
        stringify_list(STR 
            JOIN "  " 
            LINE_LENGTH 60
            HIGHLIGHT ${${VAR}}
        )
        # Wrap lines
        stringify_list(STR 
            JOIN "  * ${VAR}:\t" "\n  ...       \t"
        )
        list(APPEND ${SDO_LOG} "${STR}")
        set(${SDO_LOG} ${${SDO_LOG}} PARENT_SCOPE)
    endif()

    # Count number of selected options
    set(SELECTED_COUNT 0)
    foreach(ITEM ${${VAR}})
        list(FIND SDO_OPTIONS "${ITEM}" FOUND)
        if(${FOUND} EQUAL -1)
            # A bad option was selected.
            string(REPLACE ";" ", " TEMP "${SDO_OPTIONS}")
            dump_log(${SDO_LOG})
            message(FATAL_ERROR "\"${ITEM}\" is not valid for ${VAR} (select from: ${TEMP})")
        else()
            math(EXPR SELECTED_COUNT "${SELECTED_COUNT} + 1")
        endif()
    endforeach()

    # Check AT_LEAST rule
    if(DEFINED SDO_SELECT_AT_LEAST)
        if(${SDO_SELECT_AT_LEAST} GREATER ${SELECTED_COUNT})
            dump_log(${SDO_LOG})
            message(FATAL_ERROR "At least ${SDO_SELECT_AT_LEAST} items must be selected for ${VAR}")
        endif()
    endif()

    # Check AT_MOST rule
    if(DEFINED SDO_SELECT_AT_MOST)
        if(${SELECTED_COUNT} GREATER ${SDO_SELECT_AT_MOST})
            dump_log(${SDO_LOG})
            message(FATAL_ERROR "At most ${SDO_SELECT_AT_MOST} items can be selected for ${VAR}")
        endif()
    endif()

    # Check EXACTLY rule
    if(DEFINED SDO_SELECT_EXACTLY)
        if(NOT ${SELECTED_COUNT} EQUAL ${SDO_SELECT_EXACTLY})
            dump_log(${SDO_LOG})
            message(FATAL_ERROR "Exactly ${SDO_SELECT_EXACTLY} items must be selected for ${VAR}")
        endif()
    endif()

    # Export to parent scope
    set(${VAR} ${${VAR}} PARENT_SCOPE)
endfunction()

#[[ get_repo_version

Variable with name ${VARNAME} gets set to first 7 characters of the hash
of the last commit to the repo at ${DIR}.

Usage:
    get_repo_version(VARNAME DIR)
    
]]
macro(get_repo_version VARNAME DIR)
    execute_process(
        COMMAND git describe --tags --dirty=.dirty
        WORKING_DIRECTORY ${DIR}
        OUTPUT_VARIABLE ${VARNAME}
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
endmacro()

#[[ get_warning_suppression_flags

Returns the appropriate flags for suppressing a list of warnings, depending
on the compiler. 

Usage:
    get_warning_suppression_flags(VAR 
        [LANGUAGE <language>]
        [INTEL <warning code> ...]
        [GNU   <warning code> ...]
    )

Options:
    LANGUAGE        Returns flags will only be on for this language.

    INTEL           A list of warning codes for Intel compilers.

    GNU             A list of warnings for GCC compilers.
    
]]
function(get_warning_suppression_flags VARNAME)
    cmake_parse_arguments(DCW
        ""
        "LANGUAGE"
        "INTEL;GNU"
        ${ARGN}
    )

    if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
        if(DCW_INTEL)
            string(REPLACE ";" "," INTEL_FLAGS "${DCW_INTEL}")
            set(FLAGS "-diag-disable ${INTEL_FLAGS}")
        else()
            set(FLAGS "")
        endif()
    elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
        if(DCW_GNU)
            string(REPLACE ";" " -Wno-" GNU_FLAGS "${DCW_GNU}")
            set(FLAGS "-Wno-${GNU_FLAGS}")
        else()
            set(FLAGS "")
        endif()
    else()
        message(WARNING "Unknown compiler ID.")
        return()
    endif()

    if(DEFINED DCW_LANGUAGE)
        set(${VARNAME} "$<$<COMPILE_LANGUAGE:${DCW_LANGUAGE}>:${FLAGS}>" PARENT_SCOPE)    
    else()
        set(${VARNAME} "${FLAGS}" PARENT_SCOPE)     
    endif()
endfunction()

#[[ file_glob_directories

Returns a list of directories that match globbing expressions.

Usage:
    file_glob_directories(VAR 
        [PATTERNS <globbing-expressions> ...]
        [PATHS    <path 1> ...]
    )

Options:
    PATTERNS        Globbing expression to match. Appended to each of PATHS.

    PATHS           Paths to search. Globbing expressions are appended to each.
    
]]
function(file_glob_directories VAR)
	cmake_parse_arguments(ARGS
		"" 
		""
		"PATTERNS;PATHS" 
		${ARGN}
	)
	set(MATCHED_LIST "")
	foreach(PREFIX ${ARGS_PATHS})
        foreach(PATTERN ${ARGS_PATTERNS})
            if(IS_ABSOLUTE ${PREFIX})
                file(GLOB MATCHED ${PREFIX}/${PATTERN})
            else()
                file(GLOB MATCHED ${CMAKE_BINARY_DIR}/${PREFIX}/${PATTERN})
            endif()
            foreach(MATCHED_FILE ${MATCHED})
                get_filename_component(MATCHED_DIR ${MATCHED_FILE} DIRECTORY)
                list(APPEND MATCHED_LIST ${MATCHED_DIR})
            endforeach()
		endforeach()
    endforeach()
    if("${MATCHED_LIST}")
        list(REMOVE_DUPLICATES MATCHED_LIST)
    endif()
	set(${VAR} ${MATCHED_LIST} PARENT_SCOPE)
endfunction()

function(gc_message)
    cmake_parse_arguments(ARGS
        "" 
        "SECTION;VARIABLE;OPTION;CONSTANT;DISPLAY_NAME"
        "" 
        ${ARGN}
    )
    if(DEFINED ARGS_SECTION)
        message(STATUS "${ARGS_SECTION}:")
    endif()
    if(DEFINED ARGS_VARIABLE)
        message("  + ${ARGS_VARIABLE}: ${${ARGS_VARIABLE}}")
    endif()
    if(DEFINED ARGS_CONSTANT)
        message("    ${ARGS_DISPLAY_NAME}: ${${ARGS_CONSTANT}}")
    endif()
endfunction()