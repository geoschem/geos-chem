
#[[ add_static_target

A boilerplate function for generating static libraries in GC.

Usage:
    add_static_target(<target>
        [PUBLIC_LINK <library path | target> ...]
        [PUBLIC_INCLUDE <directory> ...]
        SOURCES <src 1> ...
    )
]]
function(add_static_target TARGET)
    # Parse arguments
    cmake_parse_arguments(ARG
        ""
        ""
        "SOURCES;PUBLIC_LINK;PUBLIC_INCLUDE"
        ${ARGN}
    )

    # Create the static library target
    add_library(${TARGET} STATIC 
        ${ARG_SOURCES}
    )

    # Direct compiled modules of ${TARGET} to include/
    set_target_properties(${TARGET}
        PROPERTIES
            Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include            
    )

    # Direct targets that depend on ${TARGET} to include module directory
    target_include_directories(${TARGET}
        PUBLIC
            ${CMAKE_CURRENT_BINARY_DIR}/include
    )

    # Set public link/include properties for ${TARGET} 
    if(DEFINED ARG_PUBLIC_LINK)
        target_link_libraries(${TARGET}
            PUBLIC
                ${ARG_PUBLIC_LINK}
        )
    endif()
    if(DEFINED ARG_PUBLIC_INCLUDE)
        target_include_directories(${TARGET}
            PUBLIC
                ${ARG_PUBLIC_INCLUDE}
        )
    endif()

    # Message that we have initialized the target's build
    message(STATUS "Initialized build target ${TARGET}")
endfunction()
