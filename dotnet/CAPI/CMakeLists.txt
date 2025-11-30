cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
set(TARGET_NAME CApi)
project(${TARGET_NAME})

file(GLOB_RECURSE SRC_FILES
    LIST_DIRECTORIES false
    *.cpp
    *.h
)

foreach(FILE_PATH IN LISTS SRC_FILES)
    file(RELATIVE_PATH REL_PATH "${CMAKE_CURRENT_SOURCE_DIR}" "${FILE_PATH}")
    string(REPLACE "/" "\\" GROUP_PATH "${REL_PATH}")
    get_filename_component(DIR_PART "${REL_PATH}" PATH)
    
    if(DIR_PART STREQUAL "")
        source_group("" FILES "${FILE_PATH}")
    else()
        string(REPLACE "/" "\\" GROUP_NAME "${DIR_PART}")
        source_group("${GROUP_NAME}" FILES "${FILE_PATH}")
    endif()
endforeach()

add_library(${TARGET_NAME} SHARED ${SRC_FILES})

target_include_directories(${TARGET_NAME} 
    PUBLIC 
        ${CMAKE_CURRENT_SOURCE_DIR}/../../src/LNLib/include
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(${TARGET_NAME} LNLib)
add_dependencies(${TARGET_NAME} LNLib)
target_compile_definitions(${TARGET_NAME} PUBLIC LNLIB_HOME)