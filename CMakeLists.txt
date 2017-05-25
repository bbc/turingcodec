cmake_minimum_required(VERSION 2.8)

if( NOT CMAKE_BUILD_TYPE )
	set( CMAKE_BUILD_TYPE Release FORCE )
endif()

project (turing)
if (USE_SYSTEM_BOOST)
else()
    add_subdirectory (boost)
endif()
add_subdirectory (havoc)
add_subdirectory (turing)
