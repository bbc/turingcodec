cmake_minimum_required(VERSION 2.8)

project(havoc)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/../../xbyak")

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
endif (MSVC)

if (CMAKE_COMPILER_IS_GNUCC)
	set(GCC_COMPILE_FLAGS "-fno-operator-names -std=c++0x -O3")
	set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COMPILE_FLAGS}")
endif (CMAKE_COMPILER_IS_GNUCC	)

add_library (havoc
	diff.h
	hadamard.h
	havoc.h
	havoc_test.h
	Jit.h
	pred_inter.h
	pred_intra.h
	quantize.h
	residual_decode.h
	sad.h
	ssd.h
	havoc.cpp
	havoc_test.c 
	diff.cpp 
	hadamard.cpp 
	pred_inter.cpp 
	pred_intra.cpp 
	quantize.cpp 
	residual_decode.cpp 
	sad.cpp 
	ssd.cpp
	README.md
)

add_executable(havoc-exe
	main.cpp
)

target_link_libraries(havoc-exe havoc)

install(TARGETS havoc LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
