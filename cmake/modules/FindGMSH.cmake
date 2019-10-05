# Module that checks whether GMSH is available.
#
# Variables used by this module which you may want to set:
# GMSH_ROOT   Path list to search for GMSH
#
# Sets the following variables
#
# GMSH_FOUND         True if GMSH was found and usable
# GMSH_INCLUDE_DIRS  Path to the GMSH include dirs
# GMSH_LIBRARIES     Name of the GMSH libraries
#

find_package(LAPACK QUIET)

#look for Gmsh.h header at positions given by the user
find_path(GMSH_INCLUDE_DIR
  NAMES "Gmsh.h"
  PATHS ${GMSH_ROOT}
  PATH_SUFFIXES "include/gmsh"
  NO_DEFAULT_PATH
)

#look for Gmsh library at positions given by the user
find_library(GMSH_LIBRARY
  NAMES "Gmsh"
  PATHS ${GMSH_ROOT}
  PATH_SUFFIXES "lib"
  NO_DEFAULT_PATH
)

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "GMSH"
  DEFAULT_MSG
  GMSH_INCLUDE_DIR
  GMSH_LIBRARY
)

mark_as_advanced(GMSH_INCLUDE_DIR GMSH_LIBRARY)

# if both header and library are found, store result
if(GMSH_FOUND)
  set(GMSH_INCLUDE_DIRS ${GMSH_INCLUDE_DIR})
  set(GMSH_LIBRARIES ${GMSH_LIBRARY})
  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determining location of GMSH succeded:\n"
    "Include directory: ${GMSH_INCLUDE_DIRS}\n"
    "Library directory: ${GMSH_LIBRARIES}\n\n")
  set(GMSH_COMPILER_FLAGS "-I${GMSH_INCLUDE_DIRS}/")
  if(LAPACK_FOUND)
    list(APPEND GMSH_LIBRARIES ${LAPACK_LIBRARIES})
  endif()
  set(GMSH_DUNE_COMPILE_FLAGS ${GMSH_COMPILER_FLAGS}
    CACHE STRING "Compile Flags used by DUNE when compiling with GMSH programs")
  set(GMSH_DUNE_LIBRARIES ${GMSH_LIBRARIES}
    CACHE STRING "Libraries used by DUNE when linking GMSH programs")
else()
  # log errornous result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKES_FILES_DIRECTORY}/CMakeError.log
    "Determing location of GMSH failed:\n"
    "Include directory: ${GMSH_INCLUDE_DIRS}\n"
    "Library directory: ${GMSH_LIBRARIES}\n\n")
endif()

#set HAVE_GMSH for config.h
set(HAVE_GMSH ${GMSH_FOUND})

# register all Gmsh related flags
if(GMSH_FOUND)
  dune_register_package_flags(
    LIBRARIES "${GMSH_LIBRARIES}"
    INCLUDE_DIRS "${GMSH_INCLUDE_DIRS}")
endif()
