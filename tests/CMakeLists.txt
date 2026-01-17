set(TARGET_NAME Tests)
project(${TARGET_NAME})
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/$<CONFIG>)
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
file(GLOB TEST_FILES ${SOURCE_DIR}/*.cpp)
add_executable(${TARGET_NAME} ${TEST_FILES})

option(INSTALL_GTEST OFF)

include(FetchContent)
FetchContent_Declare(
	googletest
	GIT_REPOSITORY https://github.com/google/googletest.git
	GIT_TAG release-1.12.0
)
FetchContent_GetProperties(googletest)
if(NOT googletest_POPULATED)
	set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
	FetchContent_MakeAvailable(googletest)	
	if(MSVC)
		foreach(_tgt gtest gtest_main gmock gmock_main)
		target_compile_definitions(${_tgt}
			PRIVATE
				"_SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING"
		)
		endforeach()
	endif()
endif()

set_target_properties(
    gtest gtest_main gmock gmock_main
    PROPERTIES FOLDER "Kernel/gtest"
)

target_include_directories(${TARGET_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/src/LNLib/include)
target_link_libraries(${TARGET_NAME} LNLib gtest gtest_main)
add_dependencies(${TARGET_NAME} LNLib)