cmake_minimum_required(VERSION 3.10)
PROJECT(MAIN)
SET(CMAKE_SHARED_LINKER_FLAGS  -lz -lm -lbz2 -llzma)

set(CMAKE_CXX_STANDARD 11)
include_directories(.)
include_directories(htslib)
link_libraries(pthread curl)
#add_custom_command(
#        OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/htslib/libhts.a
#        COMMAND make
#        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/htslib/
#        )

add_executable(linkage  src/get_relivent_snps.cpp )
add_custom_target(S ${CMAKE_CURRENT_SOURCE_DIR}/htslib/libhts.a
        COMMAND make
        WORKING-DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/htslib/
        )
target_link_libraries(linkage ${CMAKE_CURRENT_SOURCE_DIR}/htslib/libhts.a ${CMAKE_SHARED_LINKER_FLAGS})
#add_dependencies(linkage htslib/libhts)


