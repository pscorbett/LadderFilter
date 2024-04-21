cmake_minimum_required(VERSION 3.0)

set(C74_MIN_API_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../min-api)
include(${C74_MIN_API_DIR}/script/min-pretarget.cmake)


#############################################################
# MAX EXTERNAL
#############################################################

include_directories(
    "${C74_INCLUDES}"
)

set(SOURCE_FILES
    ${PROJECT_NAME}.cpp
    ../shared/signal_routing_objects.h
    ../shared/signal_routing_objects.cpp
)

add_library(${PROJECT_NAME} MODULE ${SOURCE_FILES})

# Set C++ standard to C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# Add compiler flag to disable aligned allocation
target_compile_options(${PROJECT_NAME} PRIVATE -fno-aligned-allocation)

# Set macOS deployment target for the main target
set(CMAKE_OSX_DEPLOYMENT_TARGET "14.2" CACHE STRING "Minimum macOS deployment version")

include_directories(${CMAKE_SOURCE_DIR}/source/eigen)

include(${C74_MIN_API_DIR}/script/min-posttarget.cmake)

#############################################################
# UNIT TEST
#############################################################

#include(${C74_MIN_API_DIR}/test/min-object-unittest.cmake)

