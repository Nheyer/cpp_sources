cmake_minimum_required(VERSION 3.10)
project(Recearch)


set(CMAKE_CXX_STANDARD 11)
include_directories(..)
link_libraries(pthread)
add_subdirectory(htslib)


add_executable(Recearch src/get_relivent_snps.cpp)


