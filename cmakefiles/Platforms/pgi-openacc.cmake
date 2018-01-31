### PGI Compiler for OpenMPI and OpenACC
set(TARGET_SUFFIX               ".pgi_acc")

set(ARCH                        "")
set(SIMD_SET                    "")
set(OPENMP_FLAGS                "-mp")
set(LAPACK_FLAGS                "-llapack -lblas")
set(ScaLAPACK_FLAGS             "-lscalapack -llapack -lblas")
set(ADDITIONAL_MACRO            "")
set(ADDITIONAL_OPTIMIZE_FLAGS   "")

set(Fortran_FLAGS_General       "-Mpreprocess -acc -ta=tesla,cc35,cc60,ptxinfo,maxregcount:128 -Mcuda -Minfo=acc")
set(C_FLAGS_General             "")

set(CMAKE_Fortran_COMPILER      "mpif90")
set(CMAKE_Fortran_FLAGS_DEBUG   "-pg")
set(CMAKE_Fortran_FLAGS_RELEASE "-fastsse")
set(CMAKE_C_COMPILER            "mpicc")
set(CMAKE_C_FLAGS_DEBUG         "-pg")
set(CMAKE_C_FLAGS_RELEASE       "-fastsse")

set(ENABLE_OPENACC      ON)
set(USE_MPI             ON)
# set(EXPLICIT_VEC        ON)
# set(REDUCE_FOR_MANYCORE ON)

set(LARGE_BLOCKING ON)

########
# CMake Platform-specific variables
########
set(CMAKE_SYSTEM_NAME "Linux" CACHE STRING "Cross-compiling for NVIDIA with OpenACC")
set(CMAKE_SYSTEM_PROCESSOR "NVIDIA")
