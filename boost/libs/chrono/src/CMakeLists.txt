cmake_minimum_required(VERSION 2.8)

project(boost-chrono)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../../..")

if (MSVC)
	#We statically link to reduce dependancies
	foreach(flag_var CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
		if(${flag_var} MATCHES "/MD")
			string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
		endif(${flag_var} MATCHES "/MD")
	if(${flag_var} MATCHES "/MDd")
		string(REGEX REPLACE "/MDd" "/MTd" ${flag_var} "${${flag_var}}")
	endif(${flag_var} MATCHES "/MDd")
	endforeach(flag_var)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_ALL_NO_LIB")
endif (MSVC)

if (CMAKE_COMPILER_IS_GNUCC)
	set(GCC_COMPILE_FLAGS "-fno-operator-names -std=c++0x -O3")
	set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COMPILE_FLAGS}")
endif (CMAKE_COMPILER_IS_GNUCC	)

add_library (boost_chrono
    process_cpu_clocks.cpp
	thread_clock.cpp
    chrono.cpp)
