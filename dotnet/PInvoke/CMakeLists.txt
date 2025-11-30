cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
set(TARGET_NAME LNLibSharp)
project(${TARGET_NAME})

# 检查 dotnet 是否可用
find_program(DOTNET_EXECUTABLE dotnet)
if(NOT DOTNET_EXECUTABLE)
    message(FATAL_ERROR "dotnet CLI not found. Please install .NET SDK.")
endif()

# 根据 CMake 系统名称自动推断 .NET RID（仅支持主流本地平台）
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(DOTNET_RID "linux-x64")
    # 可扩展：检测 ARM64
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64")
        set(DOTNET_RID "linux-arm64")
    endif()
elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(DOTNET_RID "win-x64")
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "ARM64")
        set(DOTNET_RID "win-arm64")
    endif()
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin") # macOS
    set(DOTNET_RID "osx-x64")
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64")
        set(DOTNET_RID "osx-arm64")
    endif()
else()
    message(WARNING "Unknown system: ${CMAKE_SYSTEM_NAME}. Defaulting to linux-x64.")
    set(DOTNET_RID "linux-x64")
endif()
message(STATUS "Detected .NET RID: ${DOTNET_RID}")

add_custom_target(${TARGET_NAME}
    COMMAND ${DOTNET_EXECUTABLE} publish
        ${CMAKE_CURRENT_SOURCE_DIR}/${TARGET_NAME}.csproj
        -r ${DOTNET_RID}
        --self-contained
        -p:PublishAot=true
        -o "${CMAKE_BINARY_DIR}/dotnet_output"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)