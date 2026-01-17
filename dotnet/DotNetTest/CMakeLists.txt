cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
set(TARGET_NAME DotNetTests)
project(${TARGET_NAME} LANGUAGES CSharp)
set(CMAKE_CSharp_FLAGS "/langversion:default")
file(GLOB_RECURSE SRC_FILES LIST_DIRECTORIES false *.cs)

add_executable(${TARGET_NAME} ${SRC_FILES})
set_property(TARGET ${TARGET_NAME} PROPERTY DOTNET_SDK "Microsoft.NET.Sdk")
set_property(TARGET ${TARGET_NAME} PROPERTY DOTNET_TARGET_FRAMEWORK_VERSION "net8.0")

set(CMAKE_CSharp_FLAGS "/platform:x64")
target_link_libraries(${TARGET_NAME} PRIVATE LNLibSharp)