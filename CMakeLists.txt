cmake_minimum_required(VERSION 3.10)
project(Recearch)
SET(CMAKE_SHARED_LINKER_FLAGS   -lz -lm -lbz2 -llzma)

set(CMAKE_CXX_STANDARD 11)
include_directories(.)
include_directories(htslib)
link_libraries(pthread curl z m bz2 lzma)
add_custom_command(OUTPUT htslib/libhts.a
        COMMAND make
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/htslib/
        )
add_executable(linkage  src/get_relivent_snps.cpp )
target_link_libraries(linkage ${CMAKE_SOURCE_DIR}/htslib/libhts.a ${CMAKE_SHARED_LINKER_FLAGS})
build_command( cp "cmake-build-debug/linkage ./")
