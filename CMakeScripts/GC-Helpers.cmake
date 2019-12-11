
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

function(gc_pretty_print)
    cmake_parse_arguments(ARGS
        "IS_BOOLEAN"
        "VARIABLE;SECTION"
        "OPTIONS"
        ${ARGN}
    )

    if(DEFINED ARGS_VARIABLE)
        if(ARGS_IS_BOOLEAN)
            set(LOGLINE "ON" "OFF")
            # Split list with "  "
            stringify_list(LOGLINE
                JOIN "  "
                LINE_LENGTH 60
            )
            # Wrap lines
            stringify_list(LOGLINE
                JOIN "  * ${ARGS_VARIABLE}:\t" "\n  ...       \t"
            )
            if("${${ARGS_VARIABLE}}")
                stringify_list(LOGLINE PRINT HIGHLIGHT "ON")
            else()
                stringify_list(LOGLINE PRINT HIGHLIGHT "OFF")
            endif()
        elseif(DEFINED ARGS_OPTIONS)
            set(LOGLINE ${ARGS_OPTIONS})
            # Split list with "  "
            stringify_list(LOGLINE
                JOIN "  "
                LINE_LENGTH 60
            )
            # Wrap lines
            stringify_list(LOGLINE
                JOIN "  * ${ARGS_VARIABLE}:\t" "\n  ...       \t"
            )
            stringify_list(LOGLINE PRINT HIGHLIGHT ${${ARGS_VARIABLE}})
        else()
            if(NOT DEFINED ${ARGS_VARIABLE})
                set(LOGLINE " ") # special case for empty variable
            else()
                set(LOGLINE ${${ARGS_VARIABLE}})
            endif()
            # Split list with "  "
            stringify_list(LOGLINE
                JOIN "  "
                LINE_LENGTH 60
            )
            # Wrap lines
            stringify_list(LOGLINE
                JOIN "  + ${ARGS_VARIABLE}:\t" "\n  ...       \t"
            )
            stringify_list(LOGLINE PRINT)
        endif()
    elseif(DEFINED ARGS_SECTION)
        message(STATUS "${ARGS_SECTION}:")
    endif()
endfunction()