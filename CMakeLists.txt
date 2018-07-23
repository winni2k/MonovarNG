# Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.7)

# Project's name
project(Monovar)

# Set make type and configs
if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")
SET(CMAKE_CXX_FLAGS " -pthread") 


# Set the output folder
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

set(CMAKE_CXX_STANDARD 14)  # enable C++14 standard


# The following folder will be included
include_directories("${PROJECT_SOURCE_DIR}/MonovarNG")
include_directories("${PROJECT_SOURCE_DIR}/MonovarNG/alglib")

# Finding htslib
find_library(HTSLIB hts)

# Finding boost
find_package(Boost REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR})
    message(STATUS ${Boost_LIBRARIES})
endif()


# Add executables
file( GLOB LIB_SOURCES ${PROJECT_SOURCE_DIR}/MonovarNG/*.cpp )
file( GLOB ALGLIB_SOURCES ${PROJECT_SOURCE_DIR}/MonovarNG/alglib/*.cpp )
# message(STATUS ${LIB_SOURCES})
# message(STATUS ${ALGLIB_SOURCES})
add_executable(monovar ${ALGLIB_SOURCES} ${LIB_SOURCES})

target_link_libraries(monovar ${HTSLIB} ${Boost_LIBRARIES})
