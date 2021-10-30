cmake_minimum_required(VERSION 3.16)
project(dgc)

find_package(Boost REQUIRED iostreams serialization)
find_package(ZLIB)
find_package(OpenMP)

add_subdirectory(vendor/libpqxx libpqxx)

option(OPTIMIZE_ULTRA "Enable all optimization of sserialize and build with lto" FALSE)
if (OPTIMIZE_ULTRA)
	set(SSERIALIZE_CONTIGUOUS_UBA_ONLY_ENABLED TRUE)
	set(SSERIALIZE_INLINE_IN_LTO_ENABLED TRUE)
	set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

if (CMAKE_INTERPROCEDURAL_OPTIMIZATION)
	set(SSERIALIZE_INLINE_IN_LTO_ENABLED TRUE)
endif()

set(SSERIALIZE_DISABLE_DTS2 TRUE)
set(SSERIALIZE_UBA_OPTIONAL_REFCOUNTING_ENABLED TRUE)
add_subdirectory(vendor/sserialize sserialize)

set(QUERY_SOURCES_CPP
	src/dgc-query.cpp
)

set(CREATE_SOURCES_CPP
	src/dgc-create.cpp
)

set(QUERY_LINK_LIBS
	ZLIB::ZLIB
	Boost::serialization
	Boost::iostreams
	pqxx
	sserialize
)

set(CREATE_LINK_LIBS
	Boost::serialization
	Boost::iostreams
	OpenMP::OpenMP_C
)

add_executable(dgc-query ${QUERY_SOURCES_CPP})
target_link_libraries(dgc-query ${QUERY_LINK_LIBS})
target_compile_features(dgc-query PRIVATE cxx_std_17)

option(DGC_SQLOA_GROUP_BOUNDS_FROM_FILE "Import group bounds from file in sqloa mode" OFF)
if (DGC_SQLOA_GROUP_BOUNDS_FROM_FILE)
	target_compile_definitions(dgc-query PRIVATE SQLOA_GROUP_BOUNDS_FROM_FILE)
endif()
option(DGC_GROUP_OF_NODE_FROM_DB "Compute group of node from db in each query in sqloa mode" OFF)
if (DGC_GROUP_OF_NODE_FROM_DB)
	target_compile_definitions(dgc-query PRIVATE GROUP_OF_NODE_FROM_DB)
endif()
option(DGC_GATHER_STATS "Gather more stats like execution times of sql queries" OFF)
if (DGC_GATHER_STATS)
	target_compile_definitions(dgc-query PRIVATE GATHER_STATS)
endif()
option(DGC_VERBOSE "Be more verbose" OFF)
if (DGC_VERBOSE)
	target_compile_definitions(dgc-query PRIVATE VERBOSE)
endif()
option(DGC_VERBOSE_QUERY "Print queries" OFF)
if (DGC_VERBOSE_QUERY)
	target_compile_definitions(dgc-query PRIVATE VERBOSE_QUERY)
endif()


add_executable(dgc-create ${CREATE_SOURCES_CPP})
target_link_libraries(dgc-create ${CREATE_LINK_LIBS})
target_compile_features(dgc-create PRIVATE cxx_std_17)
