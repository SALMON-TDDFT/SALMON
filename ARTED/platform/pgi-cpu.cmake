### PGI Compiler for OpenMPI
set(TARGET_SUFFIX               ".pgi_cpu")

set(ARCH                        "")
set(SIMD_SET                    "")
set(OPENMP_FLAGS                "-mp")
set(LAPACK_FLAGS                "-llapack -lblas")
set(ADDITIONAL_MACRO            "")
set(ADDITIONAL_OPTIMIZE_FLAGS   "")

set(Fortran_FLAGS_General       "-Mpreprocess")
set(C_FLAGS_General             "")

set(CMAKE_Fortran_COMPILER      "mpif90")
set(CMAKE_Fortran_FLAGS_DEBUG   "-pg")
set(CMAKE_Fortran_FLAGS_RELEASE "-fastsse")
set(CMAKE_C_COMPILER            "mpicc")
set(CMAKE_C_FLAGS_DEBUG         "-pg")
set(CMAKE_C_FLAGS_RELEASE       "-fastsse")

# set(STENCIL_WITH_C             ON)
# set(ENABLE_EXPLICIT_VEC        ON)
# set(ENABLE_REDUCE_FOR_MANYCORE ON)

########
# CMake Platform-specific variables
########
set(CMAKE_SYSTEM_NAME "Linux" CACHE STRING "Cross-compiling for x86-64 (Intel64 or AMD64)")
set(CMAKE_SYSTEM_PROCESSOR "x86-64")
