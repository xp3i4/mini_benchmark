# ===========================================================================
# Minibench
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (minibench CXX)

# ----------------------------------------------------------------------------
# Begin of dependency detection
# ----------------------------------------------------------------------------
#set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/seqan/include/seqan/version.h")
    set (CMAKE_INCLUDE_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}/seqan/include/
        ${CMAKE_INCLUDE_PATH})
    set (CMAKE_MODULE_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}/seqan/util/cmake
        ${CMAKE_MODULE_PATH})
endif ()
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# ----------------------------------------------------------------------------
# Add minibench targets
# ----------------------------------------------------------------------------

include_directories(include)
include_directories(external)
add_subdirectory(src)

