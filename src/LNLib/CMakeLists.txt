set(TARGET_NAME LNLib)
project(${TARGET_NAME} LANGUAGES CXX)
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
add_library(${TARGET_NAME} SHARED "")
target_compile_definitions(LNLib PRIVATE LNLIB_HOME)

target_include_directories(${TARGET_NAME} PRIVATE ${SOURCE_DIR}/include)

include(FetchContent)
FetchContent_Declare(
    eigen3
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.bz2
)
FetchContent_GetProperties(eigen3)
if (NOT eigen3_POPULATED)
  FetchContent_Populate(eigen3)
endif()

target_include_directories(${TARGET_NAME} PRIVATE ${eigen3_SOURCE_DIR})

file(GLOB rootfiles *.cpp *.h)
source_group("" FILES ${rootfiles})
target_sources(${TARGET_NAME} PRIVATE ${rootfiles})
SUBDIRLIST(SUBDIRS ${SOURCE_DIR})
foreach(subdir ${SUBDIRS})
    file(GLOB subdirFiles ${subdir}/*.cpp ${subdir}/*.h)
    string(REPLACE "/" "\\" subdir ${subdir})
    source_group(${subdir} FILES ${subdirFiles})
    target_sources(${TARGET_NAME} PRIVATE ${subdirFiles})
endforeach()