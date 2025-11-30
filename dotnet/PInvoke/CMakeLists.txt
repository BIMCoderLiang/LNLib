cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
set(TARGET_NAME LNLibSharp)
project(${TARGET_NAME} NONE)

find_program(DOTNET_EXECUTABLE dotnet)
if(NOT DOTNET_EXECUTABLE)
    message(FATAL_ERROR ".NET 8 SDK not found. Please install it.")
endif()

set(CS_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/*.cs")

if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(DOTNET_RID "win-x64")
    set(CAPI_LIB_NAME "CAPI.dll")
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "ARM64")
        set(DOTNET_RID "win-arm64")
        set(CAPI_LIB_NAME "CAPI.dll")
    endif()
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(DOTNET_RID "osx-x64")
    set(CAPI_LIB_NAME "libCAPI.dylib")
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64")
        set(DOTNET_RID "osx-arm64")
    endif()
else()
    set(DOTNET_RID "linux-x64")
    set(CAPI_LIB_NAME "libCAPI.so")
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64")
        set(DOTNET_RID "linux-arm64")
    endif()
endif()

# set(CAPI_LIB_PATH "${CMAKE_CURRENT_SOURCE_DIR}/extern/${DOTNET_RID}/${CAPI_LIB_NAME}")
# if(NOT EXISTS "${CAPI_LIB_PATH}")
#     message(FATAL_ERROR "CAPI library not found at: ${CAPI_LIB_PATH}")
# endif()

set(BUILD_DIR "${CMAKE_BINARY_DIR}")
set(CSPROJ_FILE "${BUILD_DIR}/LNLibSharp.csproj")
file(MAKE_DIRECTORY "${BUILD_DIR}")

file(WRITE "${CSPROJ_FILE}"
"<?xml version=\"1.0\" encoding=\"utf-8\"?>
<Project Sdk=\"Microsoft.NET.Sdk\">
  <PropertyGroup>
    <TargetFramework>net8.0</TargetFramework>
    <OutputType>Library</OutputType>
    <PublishAot>true</PublishAot>
    <EnableDefaultCompileItems>false</EnableDefaultCompileItems>
    <AssemblyName>LNLibSharp</AssemblyName>
    <OptimizationPreference>Size</OptimizationPreference>
  </PropertyGroup>
  <ItemGroup>
    <Compile Include=\"${CS_SOURCE}\" />
  </ItemGroup>
</Project>
")

message(STATUS "Building .NET 8 AOT plugin for RID: ${DOTNET_RID}")

set(AOT_OUTPUT_DIR "${CMAKE_BINARY_DIR}")

add_custom_target(${TARGET_NAME}
    COMMAND ${DOTNET_EXECUTABLE} publish
        "${CSPROJ_FILE}"
        -r ${DOTNET_RID}
        --self-contained
        -o "${AOT_OUTPUT_DIR}"
    DEPENDS ${CS_SOURCE}
    WORKING_DIRECTORY "${BUILD_DIR}"
    COMMENT "Publishing .NET 8 AOT plugin (calls CAPI via P/Invoke)..."
)

# add_custom_command(TARGET dotnet_aot_plugin POST_BUILD
#     COMMAND ${CMAKE_COMMAND} -E copy
#         "${CAPI_LIB_PATH}"
#         "${AOT_OUTPUT_DIR}/"
# )

# install(
#     DIRECTORY "${AOT_OUTPUT_DIR}/"
#     DESTINATION bin
# )
