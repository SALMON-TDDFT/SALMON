# recommended version : 3.0.2
cmake_minimum_required(VERSION 2.8)
enable_testing()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmakefiles/Platforms)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmakefiles/Modules)
include(${CMAKE_SOURCE_DIR}/cmakefiles/misc.cmake)
include(${CMAKE_SOURCE_DIR}/cmakefiles/setup.cmake)
include(${CMAKE_SOURCE_DIR}/cmakefiles/get_git_info.cmake)

### Project settings
project(SALMON Fortran C)

# Bulid options
## Library
option_set(USE_MPI                 "Use MPI parallelization" OFF)
option_set(USE_SCALAPACK           "Use ScaLAPACK Library"   OFF)

## Optimization
option_set(REDUCE_FOR_MANYCORE     "Enable optimization of reduction for many-core processor" OFF)
option_set(CURRENT_PREPROCESSING   "Enable preprocessing of the current computation in RT"    OFF)
option_set(LARGE_BLOCKING          "Enable large blocking for accelerator"                    OFF)

## Optimization for stencil compitations
option_set(OPT_STENCIL             "Enable optimization of the stencil"                                                ON)
option_set(DOMAIN_IS_POW2          "Enable whether the opimization assumes that a stencil domain size is power of two" OFF)
option_set(ARRAY_PADDING           "Enable array padding for the stencil"                                              ON)
option_set(EXPLICIT_VEC            "Enable explicit vectorization for the stencil"                                     OFF)
option_set(LOOP_BLOCKING           "Enable loop blocking at the stencil"                                               OFF)
option_set(SW_PREFETCH             "Enable software prefetch in the explicit vec of the stencil"                       OFF)

## Debug
option_set(USE_NVTX                "Use NVIDIA Tools Extension Library" OFF)

## test: joint hpsi subroutine
option_set(HPSI_TEST               "Use joint hpsi subroutine (test)"   OFF)
set(HPSI_TEST_SUFFIX "")
if (HPSI_TEST)
  set(HPSI_TEST_SUFFIX "_test")
endif()

## Additional library options
option_set(USE_LIBXC                "Use LibXC library"                 OFF)

# cross-compile or native-compile
if (${CMAKE_CROSSCOMPILING})
  #
  # (For supercomputer systems) cross compile mode, use toolchain file.
  #
  message(STATUS "Target platform ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
else ()
  #
  # (For PC clusters and PC) compiler environment detects automatically.
  #
  message(STATUS "Target platform is Native (${CMAKE_HOST_SYSTEM})")

  set(TARGET_SUFFIX ".cpu")

  check_mpi_compiler(${CMAKE_Fortran_COMPILER} IS_MPI_COMPILER)
  if (${IS_MPI_COMPILER})
    set(MPI_Fortran_COMPILER ${CMAKE_Fortran_COMPILER})
  endif ()

  check_mpi_compiler(${CMAKE_C_COMPILER} IS_MPI_COMPILER)
  if (${IS_MPI_COMPILER})
    set(MPI_C_COMPILER ${CMAKE_C_COMPILER})
  endif ()

  if (USE_MPI)
    find_package(MPI REQUIRED)

    if (NOT DEFINED MPI_Fortran_FOUND)
      message(FATAL_ERROR "MPI Fortran compilers not found.")
    endif()

    if (NOT DEFINED MPI_C_FOUND)
      message(FATAL_ERROR "MPI C compilers not found.")
    endif()

    set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
    set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
    
  endif ()


  set(CMAKE_Fortran_FLAGS_DEBUG   "-O2 -g")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(Fortran_FLAGS_General       "-cpp ${MPI_Fortran_COMPILE_FLAGS}")

  set(CMAKE_C_FLAGS_DEBUG         "-O2 -g")
  set(CMAKE_C_FLAGS_RELEASE       "-O3")
  set(C_FLAGS_General             "${MPI_C_COMPILE_FLAGS}")

  find_package(OpenMP REQUIRED)
  set(OPENMP_FLAGS ${OpenMP_C_FLAGS})

  if (USE_SCALAPACK)
    find_package(ScaLAPACK REQUIRED)
    set(ScaLAPACK_FLAGS ${ScaLAPACK_LINKER_FLAGS} ${ScaLAPACK_LIBRARIES})
  else ()
    find_package(LAPACK REQUIRED)
    set(LAPACK_FLAGS ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})
  endif ()
endif ()

message(STATUS "user-defined FFLAGS = ${CMAKE_Fortran_FLAGS}")
message(STATUS "user-defined CFLAGS = ${CMAKE_C_FLAGS}")

set(TARGET_NAME            "salmon${TARGET_SUFFIX}")
set(CMAKE_Fortran_FLAGS    "${ARCH} ${OPENMP_FLAGS} ${Fortran_FLAGS_General} ${ADDITIONAL_OPTIMIZE_FLAGS} ${CMAKE_Fortran_FLAGS}")
set(CMAKE_C_FLAGS          "${ARCH} ${OPENMP_FLAGS} ${C_FLAGS_General} ${ADDITIONAL_OPTIMIZE_FLAGS} ${CMAKE_C_FLAGS}")
set(EXTERNAL_LIBS          "")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR})

### MPI Setting
add_definitions_if(USE_MPI  -DSALMON_USE_MPI)

### NVIDIA Tools Extension Library
if (USE_NVTX)
  find_package(CUDA REQUIRED)
  add_definitions(-DARTED_USE_NVTX)
  link_directories(${CUDA_TOOLKIT_ROOT_DIR}/lib64)
  set(EXTERNAL_LIBS ${EXTERNAL_LIBS} nvToolsExt)
endif ()

### Numerical library
if (USE_SCALAPACK)
  set(FILE_MATHLIB scalapack)
  set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${ScaLAPACK_FLAGS})
else ()
  set(FILE_MATHLIB lapack)
  set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${LAPACK_FLAGS})
endif ()

add_definitions_if(CURRENT_PREPROCESSING -DARTED_CURRENT_PREPROCESSING)
add_definitions_if(REDUCE_FOR_MANYCORE   -DARTED_REDUCE_FOR_MANYCORE)
add_definitions_if(LARGE_BLOCKING        -DARTED_LBLK)

### Additional Libraries:
### LibXC 
if (USE_LIBXC)
  if (LIBXC_INSTALLDIR)
    include_directories(${LIBXC_INSTALLDIR}/include)
    link_directories(${LIBXC_INSTALLDIR}/lib)
    set(EXTERNAL_LIBS ${EXTERNAL_LIBS} xcf90 xc)
  else()
    find_package(LIBXCF90 REQUIRED)
    include_directories(${LIBXCF90_INCLUDES})
    set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${LIBXCF90_LIBRARIES})
  endif()
  add_definitions(-DSALMON_USE_LIBXC)
endif ()


if (OPT_STENCIL)
  add_definitions(-DARTED_STENCIL_OPTIMIZED)

  add_definitions_if(DOMAIN_IS_POW2  -DARTED_DOMAIN_POWER_OF_TWO)
  add_definitions_if(EXPLICIT_VEC    -DARTED_EXPLICIT_VECTORIZATION)
  add_definitions_if(SW_PREFETCH     -DARTED_ENABLE_SOFTWARE_PREFETCH)
  add_definitions_if(ARRAY_PADDING   -DARTED_STENCIL_PADDING)
  add_definitions_if(LOOP_BLOCKING   -DARTED_STENCIL_ENABLE_LOOP_BLOCKING)
else ()
  add_definitions(-DARTED_STENCIL_ORIGIN)
endif ()

add_definitions("${ADDITIONAL_MACRO}")

include(${CMAKE_SOURCE_DIR}/cmakefiles/workaround.cmake)


# Generate configure header
configure_file(
  ${CMAKE_SOURCE_DIR}/src/version.h.in
  ${CMAKE_BINARY_DIR}/version.h
)

configure_file(
  ${CMAKE_SOURCE_DIR}/src/versionf.h.in
  ${CMAKE_BINARY_DIR}/versionf.h
)

include_directories(${CMAKE_BINARY_DIR})


### Build target
set(ARTED_LIB            "arted_rtl")
set(GCEED_LIB            "gceed_rtl")

set(SALMON_PARALLEL_LIB  "parallel_rtl")
set(SALMON_IO_LIB        "io_rtl")
set(SALMON_MATH_LIB      "math_rtl")
set(SALMON_GS_LIB        "gs_rtl")
set(SALMON_RT_LIB        "rt_rtl")
set(SALMON_ATOM_LIB      "atom_rtl")
set(SALMON_XC_LIB        "xc_rtl")
set(SALMON_POISSON_LIB   "poisson_rtl")
set(SALMON_MAXWELL_LIB   "maxwell_rtl")
set(SALMON_COMMON_LIB    "common_rtl")
set(SALMON_EXT_LIB       "ext_rtl")
set(SALMON_MISC_LIB      "misc_rtl")

set(SALMON_SRC_DIR ${CMAKE_SOURCE_DIR}/src)
add_subdirectory(${SALMON_SRC_DIR})

# test directory
set(SALMON_TEST_DIR ${CMAKE_SOURCE_DIR}/testsuites)
add_subdirectory(${SALMON_TEST_DIR})

message(STATUS "Enabled macro")
get_directory_property(MACRO_LOG COMPILE_DEFINITIONS)
foreach(_macro ${MACRO_LOG})
  message(STATUS "  ${_macro}")
endforeach()
