cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
set(TARGET_NAME LNLibjs)
project(${TARGET_NAME})

if (CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
    message(STATUS "Building for Emscripten")
else()
    message(WARNING "Not building for Emscripten!")
endif()

add_compile_options(-std=c++17)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CONFIGURATION_TYPES "Debug;Release;RelWithDebInfo" CACHE STRING "" FORCE)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER "CMakePredefinedTargets")
set(CMAKE_EXECUTABLE_SUFFIX ".js")

file(GLOB_RECURSE SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/../src/LNLib/*.cpp")
message(STATUS "Found ${SRC_FILES} source files")
set(SOURCES ${SRC_FILES})

include_directories(${TARGET_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../build/_deps/eigen-src)
include_directories(${TARGET_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../src/LNLib/include)

add_executable(${TARGET_NAME} ${SOURCES})
target_compile_definitions(${TARGET_NAME} PUBLIC LNLIB_HOME)

target_link_options(${TARGET_NAME} PRIVATE 
    "SHELL:-s WASM=1" 
    "SHELL:-s EXPORTED_RUNTIME_METHODS=ccall,cwrap"
    "SHELL:-s MODULARIZE=1"
    "SHELL:-s EXPORT_ES6=1"
)

add_custom_command(TARGET ${TARGET_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "Build completed: $<TARGET_FILE:${TARGET_NAME}>"
)