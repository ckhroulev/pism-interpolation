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

function(find_lib PREFIX SPEC)
  find_package(PkgConfig REQUIRED)
  pkg_search_module(${PREFIX} REQUIRED IMPORTED_TARGET ${SPEC})
  if (${${PREFIX}_FOUND})
    message(STATUS "Found ${PREFIX}: ${${PREFIX}_PREFIX} (found version \"${${PREFIX}_VERSION}\")")
  endif()
endfunction()

macro(pism_find_petsc)
  # Add PETSc directories to CMAKE_PREFIX_PATH (the first on for the
  # "in place" build, the second for the "prefix" build):
  list(APPEND CMAKE_PREFIX_PATH "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}")
  list(APPEND CMAKE_PREFIX_PATH "$ENV{PETSC_DIR}")

  message(STATUS
    "Looking for PETSc (PETSC_DIR='$ENV{PETSC_DIR}', PETSC_ARCH='$ENV{PETSC_ARCH}')...")

  find_lib(PETSC "PETSc>=3.7.0")
endmacro()

pism_find_petsc()

find_lib(YAC "yac>=3.0")
find_lib(PISM "pism")
find_lib(PROJ "proj>=6.0")

add_executable(pism_interp main.cc)
target_link_libraries(pism_interp
  MPI::MPI_C
  MPI::MPI_CXX
  PkgConfig::PETSC
  PkgConfig::PISM
  PkgConfig::YAC)
