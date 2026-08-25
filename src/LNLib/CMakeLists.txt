set(TARGET_NAME LNLib)
project(${TARGET_NAME} VERSION 1.0.0 LANGUAGES CXX)
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

option(LNLib_BUILD_SHARED "Build LNLib as a shared library" ON)
if(LNLib_BUILD_SHARED)
    set(LNLib_LIBRARY_TYPE SHARED)
else()
    set(LNLib_LIBRARY_TYPE STATIC)
endif()
add_library(${TARGET_NAME} ${LNLib_LIBRARY_TYPE} "")
target_compile_definitions(LNLib PRIVATE LNLIB_HOME)
if(NOT LNLib_BUILD_SHARED)
    target_compile_definitions(LNLib PRIVATE LNLIB_STATIC)
endif()

target_include_directories(${TARGET_NAME} PRIVATE ${SOURCE_DIR}/include)

include(FetchContent)
FetchContent_Declare(
  Eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG master
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE)
set(BUILD_TESTING OFF)
set(EIGEN_BUILD_TESTING OFF)
set(EIGEN_MPL2_ONLY ON)
set(EIGEN_BUILD_PKGCONFIG OFF)
set(EIGEN_BUILD_DOC OFF)
# Keep Eigen's own install rules out of `cmake --install build` (headers and
# Eigen3Config.cmake would otherwise be copied into the install prefix).
set(EIGEN_INSTALL OFF)
FetchContent_MakeAvailable(Eigen)
target_include_directories(${TARGET_NAME} PRIVATE ${CMAKE_BINARY_DIR}/_deps/eigen-src)

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

include(CMakePackageConfigHelpers)

install(TARGETS ${TARGET_NAME}
        EXPORT LNLibTargets
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        INCLUDES DESTINATION include)

# Trailing slash on include/ is required: without it the directory itself is
# copied, producing <prefix>/include/include/.
install(DIRECTORY ${SOURCE_DIR}/include/ DESTINATION include)

install(EXPORT LNLibTargets
        FILE LNLibTargets.cmake
        NAMESPACE LNLib::
        DESTINATION lib/cmake/LNLib)

configure_package_config_file(
    ${SOURCE_DIR}/cmake/LNLibConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/LNLibConfig.cmake
    INSTALL_DESTINATION lib/cmake/LNLib)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/LNLibConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)

install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/LNLibConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/LNLibConfigVersion.cmake
        DESTINATION lib/cmake/LNLib)
