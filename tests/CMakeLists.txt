set(TERGET_NAME LNLibTests)
project(${TERGET_NAME})
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/$<CONFIG>)
add_executable(${TERGET_NAME} test.cpp)

include(FetchContent)
FetchContent_Declare(
	googletest
	GIT_REPOSITORY https://github.com/google/googletest.git
	GIT_TAG release-1.12.0
)
FetchContent_GetProperties(googletest)
if(NOT googletest_POPULATED)
	FetchContent_Populate(googletest)	
	set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
	set(gtest_disable_pthreads ON CACHE BOOL "" FORCE)	
	add_subdirectory(
		${googletest_SOURCE_DIR}
		${googletest_BINARY_DIR}
		)
		
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
    PROPERTIES FOLDER "gtest"
)

target_include_directories(${TERGET_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/src/LNLib/include)
target_link_libraries(${TERGET_NAME} LNLib gtest gtest_main)
add_dependencies(${TERGET_NAME} LNLib)