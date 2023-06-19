# ===========================================================================
# minibench: Mapping reads efficiently
# CMakeLists.txt.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (minibench)
message (STATUS "Configuring minibench/src")

find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED)
find_package (SeqAn REQUIRED)
#find_package (Eigen3 REQUIRED)
#find_package (Eigen3 3.3 REQUIRED NO_MODULE)


set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ${CXX11_CXX_FLAGS} ${SEQAN_CXX_FLAGS} \
 -Wno-class-memaccess -Wno-deprecated-copy -Wno-alloc-size-larger-than -Wno-implicit-fallthrough")

include_directories (${SEQAN_INCLUDE_DIRS})
add_definitions (${SEQAN_DEFINITIONS})
#set (CMAKE_BUILD_TYPE Debug)  
set (CMAKE_BUILD_TYPE Release)



add_library(minibench_shape_extend STATIC shape_extend.cpp)

add_executable (minibench minibench.cpp) 
target_link_libraries (minibench
                       minibench_shape_extend
                      ${SEQAN_LIBRARIES})

#
install (TARGETS minibench RUNTIME DESTINATION bin)

