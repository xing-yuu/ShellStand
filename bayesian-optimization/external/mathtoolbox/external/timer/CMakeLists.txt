cmake_minimum_required(VERSION 3.1)

project(timer CXX)
set(CMAKE_CXX_STANDARD 11)

file(GLOB header ${CMAKE_CURRENT_SOURCE_DIR}/include/timer.hpp)
install(FILES ${header} DESTINATION include)

add_library(timer INTERFACE)
target_sources(timer INTERFACE ${header})
target_include_directories(timer INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

option(TIMER_BUILD_TEST "Build test" OFF)
if(TIMER_BUILD_TEST)
	add_executable(timer-test test.cpp)
	target_link_libraries(timer-test timer)
endif()
