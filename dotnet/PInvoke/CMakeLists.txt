cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
set(TARGET_NAME LNLibSharp)
project(${TARGET_NAME} LANGUAGES CSharp)

set(CMAKE_CSharp_FLAGS "/langversion:default")

file(GLOB_RECURSE SRC_FILES
    LIST_DIRECTORIES false
    *.cs
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

set_property(TARGET ${TARGET_NAME} PROPERTY DOTNET_SDK "Microsoft.NET.Sdk")
set_property(TARGET ${TARGET_NAME} PROPERTY DOTNET_TARGET_FRAMEWORK_VERSION "net8.0")

set(CMAKE_CSharp_FLAGS "/platform:x64")
set_target_properties(${TARGET_NAME} PROPERTIES EXCLUDE_FROM_DEFAULT_BUILD TRUE)


