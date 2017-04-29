# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#.rst:
# FindBLACS
# ----------
#
# Find BLACS library
#
# This module finds an installed fortran library that implements the
# BLACS basic communication subprograms (see http://www.netlib.org/blacs/).
#
# This module is implemented that based on FindLAPACK.cmake with CMake version 3.7.2
#
# This module sets the following variables:
#
# ::
#
#   BLACS_FOUND - set to true if a library implementing the BLACS interface
#     is found
#   BLACS_LINKER_FLAGS - uncached list of required linker flags (excluding -l
#     and -L).
#   BLACS_LIBRARIES - uncached list of libraries (using full path name) to
#     link against to use BLACS
#   BLA_STATIC  if set on this determines what kind of linkage we do (static)
#   BLA_VENDOR  if set checks only the specified vendor, if not set checks
#      all the possibilities
#
# ## List of vendors (BLA_VENDOR) valid in this module
# Generic, Intel (MKL)

set(_lapack_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})

# Check the language being used
if( NOT (CMAKE_C_COMPILER_LOADED OR CMAKE_CXX_COMPILER_LOADED OR CMAKE_Fortran_COMPILER_LOADED) )
  if(BLACS_FIND_REQUIRED)
    message(FATAL_ERROR "FindBLACS requires Fortran, C, or C++ to be enabled.")
  else()
    message(STATUS "Looking for BLACS... - NOT found (Unsupported languages)")
    return()
  endif()
endif()

if (CMAKE_Fortran_COMPILER_LOADED)
  include(${CMAKE_ROOT}/Modules/CheckFortranFunctionExists.cmake)
else ()
  include(${CMAKE_ROOT}/Modules/CheckFunctionExists.cmake)
endif ()
include(${CMAKE_ROOT}/Modules/CMakePushCheckState.cmake)

cmake_push_check_state()
set(CMAKE_REQUIRED_QUIET ${BLACS_FIND_QUIETLY})

set(BLACS_FOUND FALSE)

# TODO: move this stuff to separate module

macro(Check_Lapack_Libraries LIBRARIES _prefix _name _flags _list _blas _threads)
  # This macro checks for the existence of the combination of fortran libraries
  # given by _list.  If the combination is found, this macro checks (using the
  # Check_Fortran_Function_Exists macro) whether can link against that library
  # combination using the name of a routine given by _name using the linker
  # flags given by _flags.  If the combination of libraries is found and passes
  # the link test, LIBRARIES is set to the list of complete library paths that
  # have been found.  Otherwise, LIBRARIES is set to FALSE.

  # N.B. _prefix is the prefix applied to the names of all cached variables that
  # are generated internally and marked advanced by this macro.

  set(_libraries_work TRUE)
  set(${LIBRARIES})
  set(_combined_name)
  if (NOT _libdir)
    if (WIN32)
      set(_libdir ENV LIB)
    elseif (APPLE)
      set(_libdir ENV DYLD_LIBRARY_PATH)
    else ()
      set(_libdir ENV LD_LIBRARY_PATH)
    endif ()
  endif ()
  foreach(_library ${_list})
    set(_combined_name ${_combined_name}_${_library})

    if(_libraries_work)
      if (BLA_STATIC)
        if (WIN32)
          set(CMAKE_FIND_LIBRARY_SUFFIXES .lib ${CMAKE_FIND_LIBRARY_SUFFIXES})
        endif ()
        if (APPLE)
          set(CMAKE_FIND_LIBRARY_SUFFIXES .lib ${CMAKE_FIND_LIBRARY_SUFFIXES})
        else ()
          set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
        endif ()
      else ()
        if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
          # for ubuntu's libblas3gf and liblapack3gf packages
          set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES} .so.3gf)
        endif ()
      endif ()
      find_library(${_prefix}_${_library}_LIBRARY
        NAMES ${_library}
        PATHS ${_libdir}
        )
      mark_as_advanced(${_prefix}_${_library}_LIBRARY)
      set(${LIBRARIES} ${${LIBRARIES}} ${${_prefix}_${_library}_LIBRARY})
      set(_libraries_work ${${_prefix}_${_library}_LIBRARY})
    endif()
  endforeach()

  if(_libraries_work)
    # Test this combination of libraries.
    if(UNIX AND BLA_STATIC)
      set(CMAKE_REQUIRED_LIBRARIES ${_flags} "-Wl,--start-group" ${${LIBRARIES}} ${_blas} "-Wl,--end-group" ${_threads})
    else()
      set(CMAKE_REQUIRED_LIBRARIES ${_flags} ${${LIBRARIES}} ${_blas} ${_threads})
    endif()
    #  message("DEBUG: CMAKE_REQUIRED_LIBRARIES = ${CMAKE_REQUIRED_LIBRARIES}")
    if (NOT CMAKE_Fortran_COMPILER_LOADED)
      check_function_exists("${_name}_" ${_prefix}${_combined_name}_WORKS)
    else ()
      check_fortran_function_exists(${_name} ${_prefix}${_combined_name}_WORKS)
    endif ()
    set(CMAKE_REQUIRED_LIBRARIES)
    mark_as_advanced(${_prefix}${_combined_name}_WORKS)
    set(_libraries_work ${${_prefix}${_combined_name}_WORKS})
    #message("DEBUG: ${LIBRARIES} = ${${LIBRARIES}}")
  endif()

  if(_libraries_work)
    set(${LIBRARIES} ${${LIBRARIES}} ${_blas} ${_threads})
  else()
    set(${LIBRARIES} FALSE)
  endif()

endmacro()


set(BLACS_LINKER_FLAGS)
set(BLACS_LIBRARIES)


if(BLACS_FIND_QUIETLY OR NOT BLACS_FIND_REQUIRED)
  find_package(LAPACK)
  find_package(MPI)
else()
  find_package(LAPACK REQUIRED)
  find_package(MPI    REQUIRED)
endif()


if (NOT CMAKE_Fortran_COMPILER_LOADED)
  set(_check_compiler C)
else ()
  set(_check_compiler Fortran)
endif ()


if(LAPACK_FOUND AND MPI_${_check_compiler}_FOUND)
  set(BLACS_LINKER_FLAGS "${LAPACK_LINKER_FLAGS} ${MPI_${_check_compiler}_LINK_FLAGS}")
  if (NOT $ENV{BLA_VENDOR} STREQUAL "")
    set(BLA_VENDOR $ENV{BLA_VENDOR})
  else ()
    if(NOT BLA_VENDOR)
      set(BLA_VENDOR "All")
    endif()
  endif ()

  # Generic BLACS library
  if (BLA_VENDOR STREQUAL "Generic" OR
      BLA_VENDOR STREQUAL "All")
    if ( NOT BLACS_LIBRARIES )
      check_lapack_libraries(
        BLACS_LIBRARIES
        BLACS
        pcheevx
        ""
        "scalapack"
        "${LAPACK_LIBRARIES};${MPI_${_check_compiler}_LIBRARIES}"
        ""
        )
    endif ()
  endif ()

  # Intel BLACS (MKL)
  if (BLA_VENDOR MATCHES "Intel" OR BLA_VENDOR STREQUAL "All")
    if (NOT WIN32)
      set(LM "-lm")
    endif ()
    if (CMAKE_C_COMPILER_LOADED OR CMAKE_CXX_COMPILER_LOADED)
      if(BLACS_FIND_QUIETLY OR NOT BLACS_FIND_REQUIRED)
        find_package(Threads)
      else()
        find_package(Threads REQUIRED)
      endif()

      set(BLACS_SEARCH_LIBS "")

      set(BLACS_mkl_SEARCH_SYMBOL "blacs_setup")
      set(_LIBRARIES BLACS_LIBRARIES)
      set(_LAPACK_MPI_LIBRARIES "${LAPACK_LIBRARIES};${MPI_${_check_compiler}_LIBRARIES}")

      # old
      list(APPEND BLACS_SEARCH_LIBS
        "mkl_blacs"
        "mkl_blacs_openmpi"
        "mkl_blacs_sgimpt"
        "mkl_blacs_intelmpi")
      # new >= 10.3
      list(APPEND BLACS_SEARCH_LIBS
        "mkl_blacs_lp64"
        "mkl_blacs_openmpi_lp64"
        "mkl_blacs_sgimpt_lp64"
        "mkl_blacs_intelmpi_lp64")

      # First try empty lapack libs
      if (NOT ${_LIBRARIES})
        check_lapack_libraries(
          ${_LIBRARIES}
          BLACS
          ${BLACS_mkl_SEARCH_SYMBOL}
          ""
          ""
          "${_LAPACK_MPI_LIBRARIES}"
          "${CMAKE_THREAD_LIBS_INIT};${LM}"
          )
      endif ()
      # Then try the search libs
      foreach (IT ${BLACS_SEARCH_LIBS})
        if (NOT ${_LIBRARIES})
          check_lapack_libraries(
            ${_LIBRARIES}
            BLACS
            ${BLACS_mkl_SEARCH_SYMBOL}
            ""
            "${IT}"
            "${_LAPACK_MPI_LIBRARIES}"
            "${CMAKE_THREAD_LIBS_INIT};${LM}"
            )
        endif ()
      endforeach ()
    endif ()
  endif()
else()
  message(STATUS "BLACS requires LAPACK and MPI")
endif()

if(BLACS_LIBRARIES)
  set(BLACS_FOUND TRUE)
else()
  set(BLACS_FOUND FALSE)
endif()

if(NOT BLACS_FIND_QUIETLY)
  if(BLACS_FOUND)
    message(STATUS "A library with BLACS API found.")
  else()
    if(BLACS_FIND_REQUIRED)
      message(FATAL_ERROR
        "A required library with BLACS API not found. Please specify library location."
        )
    else()
      message(STATUS
        "A library with BLACS API not found. Please specify library location."
        )
    endif()
  endif()
endif()

cmake_pop_check_state()
set(CMAKE_FIND_LIBRARY_SUFFIXES ${_lapack_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})
