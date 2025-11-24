cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
set(TARGET_NAME PyLNLib)
project(${TARGET_NAME})

add_compile_options(-std=c++17)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CONFIGURATION_TYPES "Debug;Release;RelWithDebInfo" CACHE STRING "" FORCE)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER "CMakePredefinedTargets")

find_package(Python REQUIRED COMPONENTS Interpreter Development)
message(STATUS "Python include dirs: ${Python_INCLUDE_DIRS}")

file(GLOB_RECURSE SRC_FILES src/*.cpp src/*.h)
add_library(${TARGET_NAME} SHARED ${SRC_FILES})

target_include_directories(${TARGET_NAME} 
    PUBLIC 
        ${CMAKE_CURRENT_SOURCE_DIR}/../src/LNLib/include
        ${Python_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/pybind11/include
)
target_link_libraries(${TARGET_NAME} LNLib)
add_dependencies(${TARGET_NAME} LNLib)
target_compile_definitions(${TARGET_NAME} PUBLIC LNLIB_HOME)
target_link_libraries(${TARGET_NAME} Python::Python)

set_target_properties(${TARGET_NAME} PROPERTIES SUFFIX ".pyd")
