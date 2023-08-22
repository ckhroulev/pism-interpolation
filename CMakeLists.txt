cmake_minimum_required (VERSION 3.9)

project (PISM_INTERP C CXX)

# Deal with build types
mark_as_advanced(CLEAR CMAKE_BUILD_TYPE)
# Set the build type to "Release" so that users get optimized code by default:
#
# This code block comes from the Kokkos project.
if (NOT CMAKE_BUILD_TYPE)
  set(DEFAULT_BUILD_TYPE "Release")
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING
    "Choose the type of build, options are: Debug, Release, RelWithDebInfo and MinSizeRel."
    FORCE)
endif()

# Require C++11 compiler support.
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# MPI
find_package (MPI REQUIRED COMPONENTS C CXX)
include_directories (${MPI_C_INCLUDE_PATH})

find_package(PkgConfig REQUIRED)

pkg_search_module(YAC REQUIRED IMPORTED_TARGET "yac>=3.0")
if (YAC_FOUND)
  message(STATUS
    "Found YAC: ${YAC_PREFIX} (found version \"${YAC_VERSION}\")")
endif()

pkg_search_module(PISM REQUIRED IMPORTED_TARGET "pism")
if (PISM_FOUND)
  message(STATUS
    "Found PISM: ${PISM_PREFIX} (found version \"${PISM_VERSION}\")")
endif()

add_executable(pism_interp main.cc)
target_link_libraries(pism_interp MPI::MPI_C MPI::MPI_CXX PkgConfig::PISM PkgConfig::YAC)
