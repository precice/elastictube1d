cmake_minimum_required (VERSION 3.5)

project(ElasticTube LANGUAGES CXX)
set (CMAKE_CXX_STANDARD 11)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
message(STATUS "Build configuration: " ${CMAKE_BUILD_TYPE})


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter")
endif()

find_package(precice REQUIRED CONFIG)

find_package(MPI REQUIRED
  COMPONENTS CXX)

include_directories(${MPI_CXX_INCLUDE_DIRS})
set(COMPILE_OPTIONS  ${COMPILE_FLAGS} ${MPI_CXX_COMPILE_OPTIONS})
set(LINK_FLAGS ${LINK_FLAGS} ${MPI_CXX_LINK_FLAGS})
add_definitions(${MPI_CXX_COMPILE_DEFINTIONS})

message(${MPI_CXX_LIBRARIES})

find_package(LAPACK REQUIRED)
set(LINK_FLAGS ${LINK_FLAGS} ${LAPACK_LINKER_FLAGS})


add_executable(FluidSolver
  "FluidSolver.cpp"
  "utilities.cpp"
  "FluidComputeSolution.cpp")

target_link_libraries(FluidSolver PRIVATE precice::precice)
target_link_libraries(FluidSolver PUBLIC ${LAPACK_LIBRARIES})
target_link_libraries(FluidSolver PUBLIC ${MPI_CXX_LIBRARIES})


