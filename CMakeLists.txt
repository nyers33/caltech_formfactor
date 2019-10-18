cmake_minimum_required(VERSION 2.6)
project(caltech_formfactor)

IF(WIN32)
	add_definitions(-DWIN32_LEAN_AND_MEAN)
	add_definitions(-DNOMINMAX)
	add_definitions(-D_USE_MATH_DEFINES)
ENDIF(WIN32)

add_library(caltech_formfactor STATIC
	claussen.h
	claussenp.cpp
	
	ff.h
	ffp.cpp
)

target_include_directories(caltech_formfactor PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
set_target_properties(caltech_formfactor PROPERTIES COMPILE_DEFINITIONS SINGLE)

add_executable(caltech_formfactor_demo
	main.cpp
)
target_link_libraries(caltech_formfactor_demo caltech_formfactor)

target_include_directories(caltech_formfactor_demo PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
set_target_properties(caltech_formfactor_demo PROPERTIES COMPILE_DEFINITIONS SINGLE)