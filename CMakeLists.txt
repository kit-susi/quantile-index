cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0015 NEW)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")
include(AppendCompilerFlags)
include(ExternalProject)

project(SURF CXX C)

INCLUDE_DIRECTORIES(${CMAKE_HOME_DIRECTORY}/include
                    ${CMAKE_HOME_DIRECTORY}/external/fastpfor/headers/
                    ${CMAKE_HOME_DIRECTORY}/external/cppzmq/
                    ${CMAKE_HOME_DIRECTORY}/external/zeromq/include/
                    ${CMAKE_BINARY_DIR}/external/sdsl-lite/include
                    ${CMAKE_BINARY_DIR}/external/sdsl-lite/external/libdivsufsort/include
                    )

LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/external/sdsl-lite/lib
				 ${CMAKE_BINARY_DIR}/external/zeromq/lib
					)

# Release flags
append_cxx_compiler_flags("-msse4.2 -std=c++14 -Wall" "GCC" CMAKE_CXX_FLAGS)
append_cxx_compiler_flags("-O3 -ffast-math -funroll-loops" "GCC" CMAKE_CXX_FLAGS)
append_cxx_compiler_flags("-msse4.2 -std=c++14 -funroll-loops -stdlib=libc++" "CLANG" CMAKE_CXX_FLAGS)

# Debug flags
append_cxx_compiler_flags("-O1 -g2 -fsanitize=address" "GCC" CMAKE_CXX_FLAGS_DEBUG)
append_cxx_compiler_flags("-O1 -g2 -fsanitize=address" "CLANG" CMAKE_CXX_FLAGS_DEBUG)


ADD_SUBDIRECTORY(external/zeromq)
SET_PROPERTY(DIRECTORY external/zeromq PROPERTY ZMQ_BUILD_TESTS FALSE)

ADD_SUBDIRECTORY(external/sdsl-lite)

ADD_LIBRARY(fastpfor_lib STATIC external/fastpfor/src/bitpacking.cpp
                                external/fastpfor/src/bitpackingaligned.cpp
                                external/fastpfor/src/bitpackingunaligned.cpp
                                external/fastpfor/src/simdbitpacking.cpp)

# # read the index configs
file(GLOB index_config_files RELATIVE ${CMAKE_HOME_DIRECTORY}/config/ "${CMAKE_HOME_DIRECTORY}/config/*.config")
foreach(f ${index_config_files})
	file(STRINGS ${CMAKE_HOME_DIRECTORY}/config/${f} config_contents)
	set(compile_defs "")
	foreach(keyvalue ${config_contents})
		string(REGEX REPLACE "^[ ]+" "" keyvalue ${keyvalue})
		string(REGEX MATCH "^[^=]+" key ${keyvalue})
		string(REPLACE "${key}=" "" value ${keyvalue})
		set(${key} "${value}")
		list(APPEND compile_defs ${key}=${value})
	endforeach(keyvalue)

	ADD_EXECUTABLE(surf_index-${NAME} src/surf_index.cpp)
	TARGET_LINK_LIBRARIES(surf_index-${NAME} sdsl divsufsort divsufsort64 pthread fastpfor_lib)
	set_property(TARGET surf_index-${NAME} PROPERTY COMPILE_DEFINITIONS IDXNAME="${NAME}" ${compile_defs})

    ADD_EXECUTABLE(surf_query-${NAME} src/surf_query.cpp)
    TARGET_LINK_LIBRARIES(surf_query-${NAME} sdsl divsufsort divsufsort64 pthread fastpfor_lib)
	set_property(TARGET surf_query-${NAME} PROPERTY COMPILE_DEFINITIONS IDXNAME="${NAME}" ${compile_defs})
endforeach(f)

ADD_EXECUTABLE(gen_patterns src/gen_patterns.cpp)
TARGET_LINK_LIBRARIES(gen_patterns sdsl divsufsort divsufsort64 pthread)

ADD_EXECUTABLE(gen_patterns_int src/gen_patterns_int.cpp)
TARGET_LINK_LIBRARIES(gen_patterns_int sdsl divsufsort divsufsort64 pthread)

ADD_EXECUTABLE(convert2surf src/convert2surf.cpp)
TARGET_LINK_LIBRARIES(convert2surf sdsl)

ADD_EXECUTABLE(tokenize_text src/tokenize_text.cpp)
TARGET_LINK_LIBRARIES(tokenize_text sdsl)


