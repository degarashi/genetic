cmake_minimum_required(VERSION 2.8)
project(genetic)

set(CXX_VERSION 17)
include(src/lubee/cxx_options.cmake)
enable_language(CXX)

if(NOT without-unittest)
	enable_testing()
	aux_source_directory(src TESTSRC)
	aux_source_directory(src/tests TESTSRC)
	add_executable(genetic_test ${TESTSRC})

	find_package(Threads REQUIRED)
	find_package(GTest REQUIRED)
	include_directories(
		${GTEST_INCLUDE_DIRS}
		src
	)
	target_link_libraries(genetic_test
		${GTEST_LIBRARIES}
		${GTEST_MAIN_LIBRARIES}
		${CMAKE_THREAD_LIBS_INIT}
	)
	add_test(
		NAME genetic_test
		COMMAND $<TARGET_FILE:genetic_test>
	)
endif()
